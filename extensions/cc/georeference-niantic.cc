/*
MIT License

Copyright (c) 2026 Niantic Labs

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

#include "georeference-niantic.h"
#include "load-spz.h"

#include <cmath>

namespace spz {
namespace {
template <class T>
bool readExact(std::istream& is, T& out) {
  return static_cast<bool>(is.read(reinterpret_cast<char*>(&out), sizeof(T)));
}

// Returns the number of bytes remaining between the current get position and the end of `is`,
// or -1 if the stream does not support positioning. Restores the original position on success.
std::streamoff remainingBytes(std::istream& is) {
  const std::streampos cur = is.tellg();
  if (cur == std::streampos(-1)) return -1;
  is.seekg(0, std::ios::end);
  const std::streampos end = is.tellg();
  is.seekg(cur);
  if (end == std::streampos(-1)) return -1;
  return end - cur;
}

constexpr double kUnitQuaternionTolerance = 1e-6;

// Well-formed UTF-8 (RFC 3629): no overlong forms, surrogates, or code points above U+10FFFF.
bool isValidUtf8(const std::string& s) {
  size_t i = 0;
  while (i < s.size()) {
    const auto b0 = static_cast<unsigned char>(s[i]);
    size_t len;
    uint32_t cp;
    if (b0 < 0x80) {
      ++i;
      continue;
    } else if ((b0 & 0xE0) == 0xC0) {
      len = 2;
      cp = b0 & 0x1F;
    } else if ((b0 & 0xF0) == 0xE0) {
      len = 3;
      cp = b0 & 0x0F;
    } else if ((b0 & 0xF8) == 0xF0) {
      len = 4;
      cp = b0 & 0x07;
    } else {
      return false;
    }
    if (s.size() - i < len) return false;
    for (size_t k = 1; k < len; ++k) {
      const auto b = static_cast<unsigned char>(s[i + k]);
      if ((b & 0xC0) != 0x80) return false;
      cp = (cp << 6) | (b & 0x3F);
    }
    static constexpr uint32_t kMinForLength[] = {0, 0, 0x80, 0x800, 0x10000};
    if (cp < kMinForLength[len] || cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF)) return false;
    i += len;
  }
  return true;
}

using CrsEncoding = SpzExtensionGeoreferenceNiantic::CrsEncoding;

// AuthorityCode: printable ASCII, no '@', exactly one ':', both parts non-empty.
// Projjson: well-formed UTF-8; not parsed. Both non-empty and within the u16 length field.
bool isValidCrs(const std::string& crs, CrsEncoding encoding) {
  if (crs.empty() || crs.size() > SpzExtensionGeoreferenceNiantic::kMaxCrsBytes) return false;
  if (encoding == CrsEncoding::Projjson) return isValidUtf8(crs);
  if (encoding != CrsEncoding::AuthorityCode) return false;
  size_t colons = 0;
  size_t colon = 0;
  for (size_t i = 0; i < crs.size(); ++i) {
    // Printable ASCII only: rules out spaces, control bytes and multi-byte sequences.
    const char c = crs[i];
    if (c < 0x21 || c > 0x7E) return false;
    // PROJ parses "CODE@epoch" as CoordinateMetadata; the epoch belongs in `epoch` only.
    if (c == '@') return false;
    if (c == ':') {
      ++colons;
      colon = i;
    }
  }
  return colons == 1 && colon != 0 && colon + 1 != crs.size();
}

// Fixed section: ext_version(1) + flags(1) + crs_length(2) + origin(24) + rotation(32) +
// scale(8) + epoch(8), followed by crs.
constexpr uint32_t kFixedPayloadBytes = 76;
}  // namespace

SpzExtensionGeoreferenceNiantic::SpzExtensionGeoreferenceNiantic()
    : SpzExtensionBase(SpzExtensionType::SPZ_NIANTIC_georeference) {}

uint32_t SpzExtensionGeoreferenceNiantic::payloadBytes() const {
  const uint32_t crsBytes = isValidCrs(crs, crsEncoding) ? static_cast<uint32_t>(crs.size()) : 0;
  return kFixedPayloadBytes + crsBytes;
}

void SpzExtensionGeoreferenceNiantic::write(std::ostream& os) const {
  const uint32_t t = static_cast<uint32_t>(extensionType);
  const uint32_t len = payloadBytes();
  const uint8_t extVersion = kExtVersion;
  // A transform with no target is meaningless: an invalid crs is omitted, so readers reject the record.
  const bool crsValid = isValidCrs(crs, crsEncoding);
  const uint16_t crsLength = crsValid ? static_cast<uint16_t>(crs.size()) : 0;
  if (!crsValid) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: crs (%llu bytes) is not valid for crs_encoding %u — "
           "omitted on write, readers will reject this record",
           static_cast<unsigned long long>(crs.size()), static_cast<unsigned>(crsEncoding));
  }
  uint8_t flags = 0;
  if (!std::isnan(epoch)) flags |= kFlagHasEpoch;
  flags |= static_cast<uint8_t>(static_cast<uint8_t>(crsEncoding) << kCrsEncodingShift) & kCrsEncodingMask;
  SpzLog("[SPZ] Writing extension: GeoreferenceNiantic (crs=%s)",
         !crsValid ? "<invalid>" : crsEncoding == CrsEncoding::Projjson ? "<PROJJSON>" : crs.c_str());
  os.write(reinterpret_cast<const char*>(&t), sizeof(t));
  os.write(reinterpret_cast<const char*>(&len), sizeof(len));
  os.write(reinterpret_cast<const char*>(&extVersion), sizeof(extVersion));
  os.write(reinterpret_cast<const char*>(&flags), sizeof(flags));
  os.write(reinterpret_cast<const char*>(&crsLength), sizeof(crsLength));
  os.write(reinterpret_cast<const char*>(origin.data()), sizeof(origin));
  os.write(reinterpret_cast<const char*>(rotation.data()), sizeof(rotation));
  os.write(reinterpret_cast<const char*>(&scale), sizeof(scale));
  os.write(reinterpret_cast<const char*>(&epoch), sizeof(epoch));
  if (crsLength > 0)
    os.write(crs.data(), static_cast<std::streamsize>(crsLength));
}

std::optional<SpzExtensionBasePtr> SpzExtensionGeoreferenceNiantic::read(std::istream& is) {
  SpzLog("[SPZ] Found extension: GeoreferenceNiantic");
  uint8_t extVersion{};
  uint8_t flags{};
  uint16_t crsLength{};
  auto rec = std::make_shared<SpzExtensionGeoreferenceNiantic>();
  if (!readExact(is, extVersion) || !readExact(is, flags) || !readExact(is, crsLength) ||
      !readExact(is, rec->origin) || !readExact(is, rec->rotation) || !readExact(is, rec->scale) ||
      !readExact(is, rec->epoch)) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: truncated payload — extension skipped");
    return std::nullopt;
  }
  if (extVersion != kExtVersion) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: unsupported ext_version %u (expected %u) — extension skipped",
           static_cast<unsigned>(extVersion), static_cast<unsigned>(kExtVersion));
    return std::nullopt;
  }
  // Reserved flag bits are ignored; only the bits this version defines are interpreted.
  const uint8_t crsEncoding = (flags & kCrsEncodingMask) >> kCrsEncodingShift;
  if (crsEncoding > static_cast<uint8_t>(CrsEncoding::Projjson)) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: unknown crs_encoding %u — extension skipped",
           static_cast<unsigned>(crsEncoding));
    return std::nullopt;
  }
  rec->crsEncoding = static_cast<CrsEncoding>(crsEncoding);
  const std::streamoff remaining = remainingBytes(is);
  if (remaining < 0 || crsLength > remaining) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: crs_length %u exceeds the remaining payload bytes — "
           "extension skipped",
           static_cast<unsigned>(crsLength));
    return std::nullopt;
  }
  rec->crs.resize(crsLength);
  if (crsLength > 0 && !is.read(rec->crs.data(), static_cast<std::streamsize>(crsLength))) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: failed to read crs — extension skipped");
    return std::nullopt;
  }
  if (!isValidCrs(rec->crs, rec->crsEncoding)) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: crs is empty or not valid for its crs_encoding — "
           "extension skipped");
    return std::nullopt;
  }
  for (double value : rec->origin) {
    if (!std::isfinite(value)) {
      SpzLog("[SPZ WARNING] GeoreferenceNiantic: non-finite origin — extension skipped");
      return std::nullopt;
    }
  }
  const auto& q = rec->rotation;
  const double norm = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  // Rejected rather than normalized: the library validates the record, it does not repair it.
  if (!(std::abs(norm - 1.0) <= kUnitQuaternionTolerance)) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: rotation quaternion norm %f is not unit — extension skipped",
           norm);
    return std::nullopt;
  }
  if (!std::isfinite(rec->scale) || rec->scale <= 0.0) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: scale %f is not finite and positive — extension skipped",
           rec->scale);
    return std::nullopt;
  }
  if (!(flags & kFlagHasEpoch)) {
    rec->epoch = std::numeric_limits<double>::quiet_NaN();
  } else if (!std::isfinite(rec->epoch)) {
    SpzLog("[SPZ WARNING] GeoreferenceNiantic: has_epoch is set but epoch is not finite — extension skipped");
    return std::nullopt;
  }
  return std::optional{std::move(rec)};
}

SpzExtensionType SpzExtensionGeoreferenceNiantic::type() {
  return SpzExtensionType::SPZ_NIANTIC_georeference;
}

SpzExtensionBase* SpzExtensionGeoreferenceNiantic::copyAsRawData() const {
  return new SpzExtensionGeoreferenceNiantic(*this);
}

// No PLY representation for this extension.
std::optional<SpzExtensionBasePtr> SpzExtensionGeoreferenceNiantic::tryReadFromPly(
    std::istream&, const std::unordered_set<std::string>&) const {
  return std::nullopt;
}

void SpzExtensionGeoreferenceNiantic::writePlyHeader(std::ostream&) const {}

void SpzExtensionGeoreferenceNiantic::writePlyData(std::ostream&) const {}

}  // namespace spz
