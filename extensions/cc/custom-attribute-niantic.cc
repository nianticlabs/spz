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


#include "custom-attribute-niantic.h"
#include "extension-frames.h"
#include "quantized-planes.h"
#include "load-spz.h"
#include "splat-types.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>

namespace spz {
namespace {
constexpr uint8_t kFloatBits = 32;
constexpr uint8_t kMaxSemantic = static_cast<uint8_t>(SpzAttributeSemantic::Spatial);
// Spatial attributes are made of xyz triples.
constexpr uint16_t kSpatialComponents = 3;
// Names are printable ASCII, so they are safe to log and to hand to Python and JS as strings.
constexpr char kFirstPrintable = 0x20;
constexpr char kLastPrintable = 0x7E;

size_t bytesPerComponent(uint8_t bits) {
  return isQuantizedBits(bits) ? quantizedLevelBytes(bits) : sizeof(float);
}

bool isPrintable(const std::string& s) {
  return std::all_of(s.begin(), s.end(),
                     [](char c) { return c >= kFirstPrintable && c <= kLastPrintable; });
}

void logProblem(const std::string& name, const char* problem) {
  SpzLog("[SPZ ERROR] Custom attribute '%s': %s", isPrintable(name) ? name.c_str() : "?", problem);
}

// The fields that describe an attribute in its record, shared by write and read plans.
struct Descriptor {
  std::string name;
  uint16_t dim = 0;
  SpzAttributeSemantic semantic = SpzAttributeSemantic::Custom;
  uint8_t bits = 0;
  std::vector<float> min;  // per component, only when quantized
  std::vector<float> max;
};

// Returns why the descriptor fields are invalid, or nullptr if they are valid.
const char* descriptorProblem(const std::string& name, uint16_t dim, SpzAttributeSemantic semantic,
                              uint8_t bits) {
  if (name.empty() || name.size() > std::numeric_limits<uint16_t>::max() || !isPrintable(name))
    return "name must be 1-65535 bytes of printable ASCII";
  if (dim == 0) return "dim must be at least 1";
  // An unknown semantic may need a conversion this reader can't apply.
  if (static_cast<uint8_t>(semantic) > kMaxSemantic) return "unknown semantic";
  if (semantic == SpzAttributeSemantic::Spatial && dim % kSpatialComponents != 0)
    return "SPATIAL attributes need dim to be a multiple of 3";
  if (bits != kFloatBits && !isQuantizedBits(bits)) return "bits must be 1-16 or 32";
  return nullptr;
}

// Byte size of the component-major data for numPoints points, or std::nullopt if it can't be
// addressed on this platform (e.g. on wasm32).
std::optional<size_t> dataSize(int32_t numPoints, uint16_t dim, uint8_t bits) {
  if (numPoints < 0) return std::nullopt;
  // numPoints < 2^31, dim < 2^16 and bytesPerComponent <= 4, so this fits in 64 bits.
  const uint64_t bytes = static_cast<uint64_t>(numPoints) * dim * bytesPerComponent(bits);
  if (bytes > std::numeric_limits<size_t>::max()) return std::nullopt;
  return static_cast<size_t>(bytes);
}

void computeRanges(const std::vector<float>& values, Descriptor* d) {
  d->min.assign(d->dim, 0.0f);
  d->max.assign(d->dim, 0.0f);
  const size_t points = values.size() / d->dim;
  if (points == 0) return;
  std::copy(values.begin(), values.begin() + d->dim, d->min.begin());
  std::copy(values.begin(), values.begin() + d->dim, d->max.begin());
  for (size_t p = 1; p < points; ++p) {
    for (size_t c = 0; c < d->dim; ++c) {
      d->min[c] = std::min(d->min[c], values[p * d->dim + c]);
      d->max[c] = std::max(d->max[c], values[p * d->dim + c]);
    }
  }
}

// Stores point-major values component-major (every point's first component, then every point's
// second, ...), which groups similar values together and compresses better; byte planes do the same
// for the bytes of multi-byte levels.
std::vector<uint8_t> encodeValues(const std::vector<float>& values, const Descriptor& d,
                                  size_t points, size_t bytes) {
  std::vector<uint8_t> out(bytes);
  if (points == 0) return out;
  if (isQuantizedBits(d.bits)) {
    const size_t componentBytes = points * quantizedLevelBytes(d.bits);
    for (size_t c = 0; c < d.dim; ++c) {
      encodeQuantizedPlanes(&values[c], d.dim, points, d.min[c], d.max[c], d.bits,
                            &out[c * componentBytes]);
    }
  } else {
    for (size_t c = 0; c < d.dim; ++c) {
      for (size_t p = 0; p < points; ++p) {
        std::memcpy(&out[(c * points + p) * sizeof(float)], &values[p * d.dim + c], sizeof(float));
      }
    }
  }
  return out;
}

// Inverse of encodeValues. Returns false if a quantized level exceeds 2^bits - 1.
bool decodeValues(const std::vector<uint8_t>& bytes, const Descriptor& d, size_t points,
                  std::vector<float>* values) {
  values->resize(points * d.dim);
  if (points == 0) return true;
  if (isQuantizedBits(d.bits)) {
    const size_t componentBytes = points * quantizedLevelBytes(d.bits);
    for (size_t c = 0; c < d.dim; ++c) {
      if (!decodeQuantizedPlanes(&bytes[c * componentBytes], points, d.min[c], d.max[c], d.bits,
                                 &(*values)[c], d.dim)) {
        return false;
      }
    }
  } else {
    for (size_t c = 0; c < d.dim; ++c) {
      for (size_t p = 0; p < points; ++p) {
        std::memcpy(&(*values)[p * d.dim + c], &bytes[(c * points + p) * sizeof(float)],
                    sizeof(float));
      }
    }
  }
  return true;
}

void writeRecord(std::ostream& os, const Descriptor& d,
                 const std::vector<SpzFrameLocation>& frames) {
  SpzLog("[SPZ] Writing extension: CustomAttributeNiantic (name=%s, dim=%u, bits=%u, %zu "
         "extension frames)",
         d.name.c_str(), static_cast<unsigned>(d.dim), static_cast<unsigned>(d.bits),
         frames.size());
  std::ostringstream payload;
  writeValue(payload, static_cast<uint16_t>(d.name.size()));
  payload.write(d.name.data(), static_cast<std::streamsize>(d.name.size()));
  writeValue(payload, d.dim);
  writeValue(payload, d.semantic);
  writeValue(payload, d.bits);
  payload.write(reinterpret_cast<const char*>(d.min.data()), d.min.size() * sizeof(float));
  payload.write(reinterpret_cast<const char*>(d.max.data()), d.max.size() * sizeof(float));
  writeFrameTable(payload, frames);
  writeExtensionRecord(os, SpzExtensionType::SPZ_NIANTIC_custom_attribute, payload.str());
}

// Builds the finished record from a descriptor and its decoded data, or returns nullptr (after
// logging why) if a quantized level is out of range.
SpzExtensionBasePtr resolveRecord(const Descriptor& d, int32_t numPoints,
                                  const std::vector<uint8_t>& data) {
  auto rec = std::make_shared<SpzExtensionCustomAttributeNiantic>();
  rec->name = d.name;
  rec->dim = d.dim;
  rec->semantic = d.semantic;
  rec->bits = d.bits;
  if (!decodeValues(data, d, static_cast<size_t>(numPoints), &rec->values)) {
    logProblem(d.name, "quantized level out of range");
    return nullptr;
  }
  return rec;
}
}  // namespace

SpzExtensionCustomAttributeNiantic::SpzExtensionCustomAttributeNiantic()
    : SpzExtensionBase(SpzExtensionType::SPZ_NIANTIC_custom_attribute) {}

SpzExtensionType SpzExtensionCustomAttributeNiantic::type() {
  return SpzExtensionType::SPZ_NIANTIC_custom_attribute;
}

SpzExtensionBase* SpzExtensionCustomAttributeNiantic::copyAsRawData() const {
  return new SpzExtensionCustomAttributeNiantic(*this);
}

uint32_t SpzExtensionCustomAttributeNiantic::payloadBytes() const { return 0; }

void SpzExtensionCustomAttributeNiantic::write(std::ostream& /*os*/) const {
  logProblem(name, "not written: its values can only be saved by saveSpz in SPZ version 4 files");
}

std::optional<SpzFrameWritePlan> SpzExtensionCustomAttributeNiantic::planFrameWrites(
    int32_t numPoints) const {
  if (const char* problem = descriptorProblem(name, dim, semantic, bits)) {
    logProblem(name, problem);
    return std::nullopt;
  }
  const std::optional<size_t> bytes = dataSize(numPoints, dim, bits);
  if (!bytes) {
    logProblem(name, "too large for this platform");
    return std::nullopt;
  }
  if (values.size() != static_cast<size_t>(numPoints) * dim) {
    logProblem(name, "values must hold numPoints * dim floats");
    return std::nullopt;
  }
  Descriptor d{name, dim, semantic, bits, {}, {}};
  if (isQuantizedBits(bits)) {
    if (!std::all_of(values.begin(), values.end(), [](float v) { return std::isfinite(v); })) {
      logProblem(name, "quantized values must be finite");
      return std::nullopt;
    }
    computeRanges(values, &d);
  }
  std::vector<uint8_t> data = encodeValues(values, d, static_cast<size_t>(numPoints), *bytes);
  return SpzFrameWritePlan{std::move(data),
                      [d](std::ostream& os, const std::vector<SpzFrameLocation>& frames) {
                        writeRecord(os, d, frames);
                      }};
}

std::string SpzExtensionCustomAttributeNiantic::uniquenessKey() const { return name; }

SpzExtensionBasePtr SpzExtensionCustomAttributeNiantic::withConvertedCoordinates(
    const CoordinateConverter& c) const {
  if (semantic != SpzAttributeSemantic::Spatial || isIdentityPositionTransform(c)) {
    return nullptr;
  }
  auto converted = std::make_shared<SpzExtensionCustomAttributeNiantic>(*this);
  for (size_t i = 0; i + kSpatialComponents <= converted->values.size(); i += kSpatialComponents) {
    convertDirection(c, &converted->values[i]);
  }
  return converted;
}

std::optional<SpzFrameReadPlan> parseCustomAttributeNiantic(std::istream& is) {
  Descriptor d;
  uint16_t nameLength = 0;
  if (!detail::readExact(is, nameLength)) {
    logProblem("", "truncated record");
    return std::nullopt;
  }
  d.name.resize(nameLength);
  if (!is.read(d.name.data(), nameLength) || !detail::readExact(is, d.dim) ||
      !detail::readExact(is, d.semantic) || !detail::readExact(is, d.bits)) {
    logProblem("", "truncated record");
    return std::nullopt;
  }
  if (const char* problem = descriptorProblem(d.name, d.dim, d.semantic, d.bits)) {
    logProblem(d.name, problem);
    return std::nullopt;
  }
  SpzLog("[SPZ] Found extension: CustomAttributeNiantic (name=%s)", d.name.c_str());
  if (isQuantizedBits(d.bits)) {
    d.min.resize(d.dim);
    d.max.resize(d.dim);
    if (!is.read(reinterpret_cast<char*>(d.min.data()), d.dim * sizeof(float)) ||
        !is.read(reinterpret_cast<char*>(d.max.data()), d.dim * sizeof(float))) {
      logProblem(d.name, "truncated record");
      return std::nullopt;
    }
    for (size_t c = 0; c < d.dim; ++c) {
      if (!std::isfinite(d.min[c]) || !std::isfinite(d.max[c]) || d.min[c] > d.max[c]) {
        logProblem(d.name, "invalid component range");
        return std::nullopt;
      }
    }
  }
  std::optional<std::vector<SpzFrameLocation>> frames = readFrameTable(is);
  if (!frames) {
    logProblem(d.name, "truncated frame table");
    return std::nullopt;
  }
  // Bytes after the frame table are ignored: later versions of the record may append fields.
  const uint16_t dim = d.dim;
  const uint8_t bits = d.bits;
  std::string key = d.name;
  return SpzFrameReadPlan{
    SpzExtensionType::SPZ_NIANTIC_custom_attribute, std::move(*frames), std::move(key),
    [dim, bits](int32_t numPoints) -> std::optional<uint64_t> {
      return dataSize(numPoints, dim, bits);
    },
    [d = std::move(d)](int32_t numPoints, std::vector<uint8_t>&& data) {
      return resolveRecord(d, numPoints, data);
    }};
}

std::optional<std::shared_ptr<SpzExtensionBase>> SpzExtensionCustomAttributeNiantic::tryReadFromPly(
    std::istream& /*in*/, const std::vector<PlyExtraElement>& /*elements*/) const {
  return std::nullopt;
}

void SpzExtensionCustomAttributeNiantic::writePlyHeader(std::ostream& /*out*/) const {
  SpzLog("[SPZ WARNING] Custom attribute '%s' can't be stored in PLY and was dropped",
         isPrintable(name) ? name.c_str() : "?");
}

void SpzExtensionCustomAttributeNiantic::writePlyData(std::ostream& /*out*/) const {}

std::shared_ptr<SpzExtensionCustomAttributeNiantic> findCustomAttribute(
  const std::vector<SpzExtensionBasePtr>& extensions, const std::string& name) {
  for (const auto& ext : extensions) {
    auto attr = std::dynamic_pointer_cast<SpzExtensionCustomAttributeNiantic>(ext);
    if (attr && attr->name == name) return attr;
  }
  return nullptr;
}

}  // namespace spz
