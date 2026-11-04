/*
MIT License

Copyright (c) 2026 Niantic Spatial

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include "extension-frames.h"
#include "load-spz.h"

#include <zstd.h>

#include <algorithm>
#include <cstring>
#include <limits>
#include <map>
#include <sstream>

namespace spz {
namespace {
// Extension frames are zstd skippable frames, which every zstd decoder steps over (RFC 8878
// §3.1.2), so readers that don't know an extension decode the last core stream unchanged.
// RFC 8878 reserves 0x184D2A50-0x184D2A5F for skippable frames; SPZ uses exactly this one.
constexpr uint32_t kSpzSkippableFrameMagic = 0x184D2A53;
constexpr size_t kSkippableFrameHeaderBytes = 2 * sizeof(uint32_t);  // magic, payload size
// Extension data is split into extension frames of at most this many uncompressed bytes.
// Measured at level 12 on real splat data: compression is within 0.3% of 256 MiB frames (zstd's
// match window is a few MiB, so larger frames find little more), while saving a 300 MiB attribute
// is 6x faster on 10 cores and each in-flight frame's worst-case output buffer stays small (one per
// worker; WASM runs one at a time). Not part of the format: readers accept any payload up to the
// u32 size field's 2^32 - 1 bytes.
constexpr size_t kMaxFrameUncompressedBytes = size_t{16} << 20;
static_assert(ZSTD_COMPRESSBOUND(kMaxFrameUncompressedBytes) <=
                std::numeric_limits<uint32_t>::max(),
              "a compressed payload must fit the skippable header's u32 size field");

using RecordKey = std::pair<uint32_t, std::string>;  // extension type, uniqueness key

uint32_t readU32(const uint8_t* p) {
  uint32_t v = 0;
  std::memcpy(&v, p, sizeof(v));
  return v;
}

void logDropped(SpzExtensionType type, const char* reason) {
  SpzLog("[SPZ ERROR] Dropped extension 0x%08x: %s", static_cast<unsigned>(type), reason);
}

std::optional<RecordKey> keyOf(const SpzExtensionBasePtr& record) {
  const auto* keyed = dynamic_cast<const SpzKeyedExtension*>(record.get());
  if (!keyed) return std::nullopt;
  return RecordKey{static_cast<uint32_t>(record->extensionType), keyed->uniquenessKey()};
}

std::optional<RecordKey> keyOf(const SpzFrameReadPlan& readPlan) {
  if (!readPlan.key) return std::nullopt;
  return RecordKey{static_cast<uint32_t>(readPlan.type), *readPlan.key};
}

// Returns one flag per key: nonzero when another entry has the same key.
std::vector<uint8_t> findDuplicates(const std::vector<std::optional<RecordKey>>& keys) {
  std::map<RecordKey, size_t> counts;
  for (const auto& key : keys) {
    if (key) counts[*key]++;
  }
  std::vector<uint8_t> duplicate(keys.size(), 0);
  for (size_t i = 0; i < keys.size(); i++) {
    duplicate[i] = keys[i] && counts[*keys[i]] > 1;
  }
  return duplicate;
}

// Returns why loc doesn't name a usable extension frame after the zstd frame of the host
// stream, or nullptr if it does.
const char* frameProblem(const SpzFrameLocation& loc, const uint8_t* host, size_t hostSize,
                         size_t coreFrameSize) {
  if (loc.frameOffset < coreFrameSize || loc.frameOffset > hostSize ||
      hostSize - loc.frameOffset < kSkippableFrameHeaderBytes ||
      loc.compressedSize > hostSize - loc.frameOffset - kSkippableFrameHeaderBytes) {
    return "frame lies outside the last core stream";
  }
  const uint8_t* frame = host + loc.frameOffset;
  if (readU32(frame) != kSpzSkippableFrameMagic ||
      readU32(frame + sizeof(uint32_t)) != loc.compressedSize) {
    return "frame header doesn't match the record";
  }
  const uint8_t* payload = frame + kSkippableFrameHeaderBytes;
  const size_t payloadSize = static_cast<size_t>(loc.compressedSize);
  if (ZSTD_findFrameCompressedSize(payload, payloadSize) != payloadSize) {
    return "frame payload isn't exactly one zstd frame";
  }
  // The zstd frame header records its decompressed size; requiring it to match before allocating
  // stops a record from claiming sizes its data doesn't have.
  const unsigned long long contentSize = ZSTD_getFrameContentSize(payload, payloadSize);
  if (contentSize == ZSTD_CONTENTSIZE_UNKNOWN || contentSize == ZSTD_CONTENTSIZE_ERROR ||
      loc.uncompressedSize == 0 || contentSize != loc.uncompressedSize) {
    return "frame size doesn't match its zstd frame";
  }
  return nullptr;
}

// Returns why readPlan's extension frames can't be decoded from host, or nullptr after setting
// *dataBytes to the record's decoded size.
const char* framesProblem(const SpzFrameReadPlan& readPlan, int32_t numPoints,
                          const uint8_t* host, size_t hostSize, size_t coreFrameSize,
                          uint64_t* dataBytes) {
  const std::optional<uint64_t> expected = readPlan.dataBytes(numPoints);
  if (!expected || *expected > std::numeric_limits<size_t>::max()) {
    return "data too large for this platform";
  }
  if (ZSTD_isError(coreFrameSize)) return "the last core stream doesn't start with a zstd frame";
  uint64_t total = 0;
  for (const SpzFrameLocation& loc : readPlan.frames) {
    if (const char* problem = frameProblem(loc, host, hostSize, coreFrameSize)) return problem;
    if (loc.uncompressedSize > *expected - total) return "frames hold more data than the record";
    total += loc.uncompressedSize;
  }
  if (total != *expected) return "frames hold less data than the record";
  *dataBytes = total;
  return nullptr;
}
}  // namespace

void writeExtensionRecord(std::ostream& os, SpzExtensionType type, const std::string& payload) {
  writeValue(os, static_cast<uint32_t>(type));
  writeValue(os, static_cast<uint32_t>(payload.size()));
  os.write(payload.data(), static_cast<std::streamsize>(payload.size()));
}

void writeFrameTable(std::ostream& os, const std::vector<SpzFrameLocation>& frames) {
  writeValue(os, static_cast<uint32_t>(frames.size()));
  for (const SpzFrameLocation& f : frames) {
    writeValue(os, f.frameOffset);
    writeValue(os, f.compressedSize);
    writeValue(os, f.uncompressedSize);
  }
}

std::optional<std::vector<SpzFrameLocation>> readFrameTable(std::istream& is) {
  uint32_t count = 0;
  if (!detail::readExact(is, count)) return std::nullopt;
  std::vector<SpzFrameLocation> frames;
  for (uint32_t i = 0; i < count; ++i) {
    SpzFrameLocation f;
    if (!detail::readExact(is, f.frameOffset) || !detail::readExact(is, f.compressedSize) ||
        !detail::readExact(is, f.uncompressedSize)) {
      return std::nullopt;
    }
    frames.push_back(f);
  }
  return frames;
}

std::optional<ExtensionFrameEncoder> ExtensionFrameEncoder::plan(
    const std::vector<SpzExtensionBasePtr>& extensions, int32_t numPoints) {
  ExtensionFrameEncoder encoder;
  for (const SpzExtensionBasePtr& ext : extensions) {
    if (!ext) continue;
    Entry entry{ext, std::nullopt, 0, 0};
    if (const auto* framed = dynamic_cast<const SpzFramedExtension*>(ext.get())) {
      entry.writePlan = framed->planFrameWrites(numPoints);
      if (!entry.writePlan) {
        SpzLog("[SPZ ERROR] saveSpz: extension 0x%08x can't be saved",
               static_cast<unsigned>(ext->extensionType));
        return std::nullopt;
      }
      // The plan's data lives on the heap, so these pointers survive moving entries and the
      // encoder.
      const std::vector<uint8_t>& data = entry.writePlan->data;
      entry.firstFrame = encoder.frameData_.size();
      for (size_t offset = 0; offset < data.size(); offset += kMaxFrameUncompressedBytes) {
        encoder.frameData_.push_back(
          {data.data() + offset, std::min(kMaxFrameUncompressedBytes, data.size() - offset)});
      }
      entry.frameCount = encoder.frameData_.size() - entry.firstFrame;
    }
    encoder.entries_.push_back(std::move(entry));
  }
  std::vector<std::optional<RecordKey>> keys;
  for (const Entry& entry : encoder.entries_) {
    keys.push_back(keyOf(entry.record));
  }
  const std::vector<uint8_t> duplicate = findDuplicates(keys);
  for (size_t i = 0; i < keys.size(); i++) {
    if (duplicate[i]) {
      SpzLog("[SPZ ERROR] saveSpz: two extensions of type 0x%08x share the key '%s'",
             static_cast<unsigned>(keys[i]->first), keys[i]->second.c_str());
      return std::nullopt;
    }
  }
  return encoder;
}

std::string ExtensionFrameEncoder::writeRecords(const std::vector<std::vector<uint8_t>>& compressed,
                                                uint64_t firstFrameOffset) const {
  std::ostringstream os;
  uint64_t frameOffset = firstFrameOffset;
  for (const Entry& entry : entries_) {
    if (!entry.writePlan) {
      entry.record->write(os);
      continue;
    }
    std::vector<SpzFrameLocation> frames;
    for (size_t k = entry.firstFrame; k < entry.firstFrame + entry.frameCount; k++) {
      frames.push_back({frameOffset, compressed[k].size(), frameData_[k].size});
      frameOffset += kSkippableFrameHeaderBytes + compressed[k].size();
    }
    entry.writePlan->writeRecord(os, frames);
  }
  return os.str();
}

uint64_t skippableFramesBytes(const std::vector<std::vector<uint8_t>>& compressed) {
  uint64_t bytes = 0;
  for (const auto& chunk : compressed) bytes += kSkippableFrameHeaderBytes + chunk.size();
  return bytes;
}

void appendSkippableFrames(const std::vector<std::vector<uint8_t>>& compressed,
                           std::vector<uint8_t>* out) {
  for (const auto& chunk : compressed) {
    const uint32_t header[2] = {kSpzSkippableFrameMagic, static_cast<uint32_t>(chunk.size())};
    const auto* headerBytes = reinterpret_cast<const uint8_t*>(header);
    out->insert(out->end(), headerBytes, headerBytes + sizeof(header));
    out->insert(out->end(), chunk.begin(), chunk.end());
  }
}

ExtensionFrameDecoder ExtensionFrameDecoder::fromZone(const uint8_t* zone, size_t zoneSize,
                                                      int32_t numPoints, const uint8_t* host,
                                                      size_t hostSize, uint64_t maxDataBytes,
                                                      const SpzFrameReadFilter& shouldDecode) {
  std::istringstream stream(std::string(reinterpret_cast<const char*>(zone), zoneSize));
  std::vector<SpzParsedExtension> parsed;
  readExtensionRecords(stream, parsed);
  return ExtensionFrameDecoder(std::move(parsed), numPoints, host, hostSize, maxDataBytes,
                               shouldDecode);
}

ExtensionFrameDecoder::ExtensionFrameDecoder(std::vector<SpzParsedExtension> parsed,
                                             int32_t numPoints, const uint8_t* host,
                                             size_t hostSize, uint64_t maxDataBytes,
                                             const SpzFrameReadFilter& shouldDecode)
    : numPoints_(numPoints) {
  // Duplicates are dropped before any extension frame is checked, so every reader agrees on which
  // records survive whatever state their extension frames are in.
  std::vector<std::optional<RecordKey>> keys;
  for (const SpzParsedExtension& p : parsed) {
    keys.push_back(p.record ? keyOf(p.record) : p.readPlan ? keyOf(*p.readPlan) : std::nullopt);
  }
  const std::vector<uint8_t> duplicate = findDuplicates(keys);

  std::optional<size_t> coreFrameSize;
  uint64_t totalDataBytes = 0;
  // Reserved so entries never move: tasks point into their data buffers.
  entries_.reserve(parsed.size());
  for (size_t i = 0; i < parsed.size(); i++) {
    SpzParsedExtension& p = parsed[i];
    if (duplicate[i]) {
      logDropped(p.record ? p.record->extensionType : p.readPlan->type,
                 "another record has the same key");
      continue;
    }
    if (p.record) {
      entries_.push_back({std::move(p.record), std::nullopt, {}, 0});
      continue;
    }
    if (!p.readPlan) continue;
    // After the duplicate check, so filtering never changes which records survive.
    if (!shouldDecode(*p.readPlan)) {
      SpzLog("[SPZ] Skipped extension 0x%08x%s%s: not requested",
             static_cast<unsigned>(p.readPlan->type), p.readPlan->key ? " " : "",
             p.readPlan->key ? p.readPlan->key->c_str() : "");
      continue;
    }
    if (!coreFrameSize) coreFrameSize = ZSTD_findFrameCompressedSize(host, hostSize);
    uint64_t dataBytes = 0;
    if (const char* problem =
          framesProblem(*p.readPlan, numPoints, host, hostSize, *coreFrameSize, &dataBytes)) {
      logDropped(p.readPlan->type, problem);
      continue;
    }
    if (dataBytes > maxDataBytes - totalDataBytes) {
      logDropped(p.readPlan->type, "more extension data than the file size allows");
      continue;
    }
    totalDataBytes += dataBytes;
    entries_.push_back({nullptr, std::move(p.readPlan), {}, tasks_.size()});
    Entry& entry = entries_.back();
    entry.data.resize(static_cast<size_t>(dataBytes));
    size_t written = 0;
    for (const SpzFrameLocation& loc : entry.readPlan->frames) {
      tasks_.push_back({host + loc.frameOffset + kSkippableFrameHeaderBytes,
                        static_cast<size_t>(loc.compressedSize), entry.data.data() + written,
                        static_cast<size_t>(loc.uncompressedSize)});
      written += static_cast<size_t>(loc.uncompressedSize);
    }
  }
}

std::vector<SpzExtensionBasePtr> ExtensionFrameDecoder::finish(
    const std::vector<uint8_t>& succeeded) && {
  std::vector<SpzExtensionBasePtr> records;
  for (Entry& entry : entries_) {
    if (entry.record) {
      records.push_back(std::move(entry.record));
      continue;
    }
    const auto first = succeeded.begin() + static_cast<std::ptrdiff_t>(entry.firstTask);
    const auto last = first + static_cast<std::ptrdiff_t>(entry.readPlan->frames.size());
    if (std::find(first, last, 0) != last) {
      logDropped(entry.readPlan->type, "frame decompression failed");
      continue;
    }
    // resolve() logs why it rejects the data.
    if (SpzExtensionBasePtr record = entry.readPlan->resolve(numPoints_, std::move(entry.data))) {
      records.push_back(std::move(record));
    }
  }
  return records;
}

void dropRecordsWithDuplicateKeys(std::vector<SpzExtensionBasePtr>& extensions) {
  extensions.erase(std::remove(extensions.begin(), extensions.end(), nullptr), extensions.end());
  std::vector<std::optional<RecordKey>> keys;
  for (const auto& ext : extensions) keys.push_back(keyOf(ext));
  const std::vector<uint8_t> duplicate = findDuplicates(keys);
  std::vector<SpzExtensionBasePtr> kept;
  for (size_t i = 0; i < extensions.size(); i++) {
    if (duplicate[i]) {
      logDropped(extensions[i]->extensionType, "another record has the same key");
    } else {
      kept.push_back(std::move(extensions[i]));
    }
  }
  extensions = std::move(kept);
}

bool keysAreUnique(const std::vector<SpzExtensionBasePtr>& extensions) {
  std::vector<std::optional<RecordKey>> keys;
  for (const auto& ext : extensions) keys.push_back(ext ? keyOf(ext) : std::nullopt);
  const std::vector<uint8_t> duplicate = findDuplicates(keys);
  for (size_t i = 0; i < keys.size(); i++) {
    if (duplicate[i]) {
      SpzLog("[SPZ ERROR] two extensions of type 0x%08x share the key '%s'",
             static_cast<unsigned>(keys[i]->first), keys[i]->second.c_str());
      return false;
    }
  }
  return true;
}

}  // namespace spz
