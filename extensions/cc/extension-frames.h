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


// Internal to the library (not installed): stores framed extensions' data (SpzFramedExtension) in
// extension frames after the last core stream of NGSP (v4) files, and reads it back.

#ifndef SPZ_EXTENSION_FRAMES_H_
#define SPZ_EXTENSION_FRAMES_H_

#include "splat-extensions.h"
#include "zstd-batch.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace spz {

// How to load a framed record: which extension frames hold its data, and how to build the record
// once the library has decoded them. Built-in parsers return these from a record's payload.
struct SpzFrameReadPlan {
  // The type of the record this resolves to.
  SpzExtensionType type;
  // The extension frames holding the record's data, in order.
  std::vector<SpzFrameLocation> frames;
  // The key the finished record's SpzKeyedExtension::uniquenessKey() will return, for types that
  // have one, so duplicates are dropped before anything is decoded.
  std::optional<std::string> key;
  // Size the decoded data must have for a cloud of numPoints points, or std::nullopt if it can't be
  // addressed on this platform.
  std::function<std::optional<uint64_t>(int32_t numPoints)> dataBytes;
  // Builds the record from its extension frames' decompressed payloads, concatenated in order.
  // Returns nullptr (after logging why) if the data is invalid. Precondition:
  // data.size() == dataBytes(numPoints).
  std::function<SpzExtensionBasePtr(int32_t numPoints, std::vector<uint8_t>&& data)> resolve;
};

// Whether the caller wants a framed record decoded.
using SpzFrameReadFilter = std::function<bool(const SpzFrameReadPlan& readPlan)>;

// One record read from an extension zone: a finished record, or a framed record still waiting for
// its extension frames. Exactly one of the two is set.
struct SpzParsedExtension {
  SpzExtensionBasePtr record;
  std::optional<SpzFrameReadPlan> readPlan;
};

// Reads every record from an extension zone (defined in splat-extensions.cc, next to the parsers).
void readExtensionRecords(std::istream& is, std::vector<SpzParsedExtension>& out);

// Built-in parsers for framed records, called by readExtensionRecords with a record's payload (the
// bytes after type and byteLength). Each returns std::nullopt (after logging why) if the payload is
// malformed.
std::optional<SpzFrameReadPlan> parseCustomAttributeNiantic(std::istream& payload);

// ---- Helpers for writing and parsing framed extensions' records ----

// Writes a trivially copyable value in the file's (little-endian, native) byte order; the reading
// counterpart is detail::readExact in splat-extensions.h.
template <class T>
void writeValue(std::ostream& os, const T& v) {
  os.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

// Writes a whole extension record: type, the payload's byte length, then the payload.
// Precondition: payload is at most 2^32 - 1 bytes.
void writeExtensionRecord(std::ostream& os, SpzExtensionType type, const std::string& payload);

// Writes a frame table: u32 count, then {u64 frameOffset, u64 compressedSize, u64 uncompressedSize}
// per extension frame.
void writeFrameTable(std::ostream& os, const std::vector<SpzFrameLocation>& frames);

// Reads a frame table written by writeFrameTable, or returns std::nullopt if it's truncated.
// Entries are read one at a time, so a count that lies can't drive a large allocation.
std::optional<std::vector<SpzFrameLocation>> readFrameTable(std::istream& is);

// Plans where every extension record's data goes when a cloud is saved.
class ExtensionFrameEncoder {
 public:
  // Plans the records for a cloud of numPoints points. Returns std::nullopt (after logging why) if
  // a framed record is invalid or two records share a uniqueness key.
  static std::optional<ExtensionFrameEncoder> plan(
    const std::vector<SpzExtensionBasePtr>& extensions, int32_t numPoints);

  bool hasFrames() const { return !frameData_.empty(); }

  // Uncompressed contents of every extension frame, in file order. Valid for the encoder's
  // lifetime.
  const std::vector<ByteSpan>& frameData() const { return frameData_; }

  // Returns the extension zone: every record, framed records listing their extension frames.
  // Preconditions: compressed holds frameData(), each compressed into one zstd frame, in order; the
  // extension frames will be appended right after the last core stream's zstd frame of
  // firstFrameOffset bytes.
  std::string writeRecords(const std::vector<std::vector<uint8_t>>& compressed,
                           uint64_t firstFrameOffset) const;

 private:
  struct Entry {
    SpzExtensionBasePtr record;
    std::optional<SpzFrameWritePlan> writePlan;  // unset for records without extension frames
    size_t firstFrame = 0;
    size_t frameCount = 0;
  };
  std::vector<Entry> entries_;
  std::vector<ByteSpan> frameData_;
};

// Bytes the compressed chunks take once wrapped as extension frames.
uint64_t skippableFramesBytes(const std::vector<std::vector<uint8_t>>& compressed);

// Appends the compressed chunks to out as extension frames. Precondition: each chunk is at most
// 2^32 - 1 bytes (guaranteed for data planned by ExtensionFrameEncoder).
void appendSkippableFrames(const std::vector<std::vector<uint8_t>>& compressed,
                           std::vector<uint8_t>* out);

// Turns parsed extension records back into finished records: validates each framed record's
// extension frames against the file, hands out the decompression work, then resolves the records.
class ExtensionFrameDecoder {
 public:
  // Parses an extension zone (zoneSize bytes) and prepares its framed records for decoding. host
  // is the last core stream (hostSize bytes); every pointer is into the loaded file, which must
  // outlive the decoder.
  static ExtensionFrameDecoder fromZone(const uint8_t* zone, size_t zoneSize, int32_t numPoints,
                                        const uint8_t* host, size_t hostSize,
                                        uint64_t maxDataBytes,
                                        const SpzFrameReadFilter& shouldDecode);

  // host is the last core stream (hostSize bytes, inside the loaded file). maxDataBytes caps the
  // total decoded size of all framed records. Records sharing a uniqueness key, and framed records
  // whose extension frames are invalid, are dropped (logged) here, before anything is allocated.
  // Framed records shouldDecode rejects are skipped unread. Precondition: shouldDecode is
  // callable.
  ExtensionFrameDecoder(std::vector<SpzParsedExtension> parsed, int32_t numPoints,
                        const uint8_t* host, size_t hostSize, uint64_t maxDataBytes,
                        const SpzFrameReadFilter& shouldDecode);

  // Decompression work for the remaining framed records. The destinations are owned by the
  // decoder, which must outlive the tasks.
  const std::vector<DecompressTask>& tasks() const { return tasks_; }

  // Returns the finished records in file order, consuming the decoder. Records whose extension
  // frames failed to decompress or whose data is invalid are dropped. Precondition: succeeded[i]
  // says whether tasks()[i] produced exactly its dstSize bytes.
  std::vector<SpzExtensionBasePtr> finish(const std::vector<uint8_t>& succeeded) &&;

 private:
  struct Entry {
    SpzExtensionBasePtr record;                // set for records without extension frames
    std::optional<SpzFrameReadPlan> readPlan;  // set for framed records
    std::vector<uint8_t> data;                     // the framed record's decoded data
    size_t firstTask = 0;
  };
  int32_t numPoints_;
  std::vector<Entry> entries_;
  std::vector<DecompressTask> tasks_;
};

// Drops every record whose uniqueness key another record of the same type shares (logged).
void dropRecordsWithDuplicateKeys(std::vector<SpzExtensionBasePtr>& extensions);

// False (after logging the key) if two records of the same type share a uniqueness key.
bool keysAreUnique(const std::vector<SpzExtensionBasePtr>& extensions);

}  // namespace spz

#endif  // SPZ_EXTENSION_FRAMES_H_
