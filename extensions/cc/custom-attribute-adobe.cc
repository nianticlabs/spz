/*
MIT License

Copyright (c) 2026 Adobe Inc.

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

#include "custom-attribute-adobe.h"
#include "load-spz.h"

#include <zstd.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace spz {
namespace {
constexpr uint8_t kFloatBits = 32;
// On 6-bit attributes of a 500k-splat asset, level 3 came out 2% smaller and 5x faster than 12.
constexpr int kZstdLevel = 3;
// A 4-byte RLE block expands to at most 128 KiB, so no zstd frame inflates more than 32768x.
// Rejecting larger claims bounds the allocation for malformed files.
constexpr uint64_t kMaxZstdExpansion = 32768;

template <class T>
bool readExact(std::istream& is, T& out) {
  return static_cast<bool>(is.read(reinterpret_cast<char*>(&out), sizeof(T)));
}

template <class T>
void append(std::string& out, const T* data, size_t count) {
  out.append(reinterpret_cast<const char*>(data), count * sizeof(T));
}

bool isQuantized(uint8_t bits) { return bits >= 1 && bits <= 8; }

// nameLength, name, dim, kind, bits, then min[dim] and max[dim] when quantized.
uint64_t settingsBytes(size_t nameLength, uint16_t dim, uint8_t bits) {
  return sizeof(uint16_t) + nameLength + sizeof(uint16_t) + 2 +
         (isQuantized(bits) ? 2 * sizeof(float) * dim : 0);
}
}  // namespace

SpzExtensionCustomAttributeAdobe::SpzExtensionCustomAttributeAdobe()
    : SpzExtensionBase(SpzExtensionType::SPZ_ADOBE_custom_attribute) {}

bool SpzExtensionCustomAttributeAdobe::validate(int32_t numPoints) const {
  const size_t rawBytes = values.size() * (isQuantized(bits) ? 1 : sizeof(float));
  const char* problem = nullptr;
  if (name.empty() || name.size() > std::numeric_limits<uint16_t>::max())
    problem = "name must be 1-65535 bytes";
  else if (dim == 0)
    problem = "dim must be at least 1";
  else if (kind != SpzCustomAttributeKind::Scalar && kind != SpzCustomAttributeKind::Vector3)
    problem = "unknown kind";
  else if (kind == SpzCustomAttributeKind::Vector3 && dim % 3 != 0)
    problem = "a Vector3 attribute needs dim to be a multiple of 3";
  else if (bits != kFloatBits && !isQuantized(bits))
    problem = "bits must be 1-8 or 32";
  else if (values.size() != static_cast<size_t>(numPoints) * dim)
    problem = "values must hold numPoints * dim floats";
  else if (isQuantized(bits) &&
           !std::all_of(values.begin(), values.end(), [](float v) { return std::isfinite(v); }))
    problem = "quantized values must be finite";
  else if (settingsBytes(name.size(), dim, bits) + ZSTD_compressBound(rawBytes) >
           std::numeric_limits<uint32_t>::max())
    problem = "data exceeds the 4 GiB record limit";
  if (problem) {
    SpzLog("[SPZ ERROR] Custom attribute '%s': %s", name.c_str(), problem);
    return false;
  }
  return true;
}

std::string SpzExtensionCustomAttributeAdobe::encode() const {
  const uint16_t nameLength = static_cast<uint16_t>(name.size());
  const uint8_t kindValue = static_cast<uint8_t>(kind);
  std::string out;
  append(out, &nameLength, 1);
  out.append(name, 0, nameLength);
  append(out, &dim, 1);
  append(out, &kindValue, 1);
  append(out, &bits, 1);

  // The file stores the values component by component (every point's first value, then every
  // point's second, ...): on real 6-bit data that is 7% smaller than point by point.
  const size_t points = dim > 0 ? values.size() / dim : 0;
  const void* src = values.data();  // dim == 1 is already in file order
  size_t srcBytes = values.size() * sizeof(float);
  std::vector<uint8_t> levels;
  std::vector<float> planar;
  if (isQuantized(bits) && dim > 0) {
    std::vector<float> lo(dim, 0.0f), hi(dim, 0.0f);
    if (points > 0) {
      std::copy(values.begin(), values.begin() + dim, lo.begin());
      std::copy(values.begin(), values.begin() + dim, hi.begin());
    }
    for (size_t p = 0; p < points; ++p) {
      for (size_t c = 0; c < dim; ++c) {
        lo[c] = std::min(lo[c], values[p * dim + c]);
        hi[c] = std::max(hi[c], values[p * dim + c]);
      }
    }
    append(out, lo.data(), dim);
    append(out, hi.data(), dim);

    // Double precision keeps the range finite even when hi - lo overflows float.
    const double maxLevel = static_cast<double>((1 << bits) - 1);
    std::vector<double> scale(dim);
    for (size_t c = 0; c < dim; ++c)
      scale[c] = hi[c] > lo[c] ? maxLevel / (static_cast<double>(hi[c]) - lo[c]) : 0.0;
    levels.resize(points * dim);
    for (size_t p = 0; p < points; ++p) {
      for (size_t c = 0; c < dim; ++c) {
        const double level = std::round((values[p * dim + c] - static_cast<double>(lo[c])) * scale[c]);
        levels[c * points + p] = static_cast<uint8_t>(std::min(level, maxLevel));
      }
    }
    src = levels.data();
    srcBytes = levels.size();
  } else if (dim > 1) {
    planar.resize(points * dim);
    for (size_t p = 0; p < points; ++p) {
      for (size_t c = 0; c < dim; ++c)
        planar[c * points + p] = values[p * dim + c];
    }
    src = planar.data();
    srcBytes = planar.size() * sizeof(float);
  }

  const size_t offset = out.size();
  out.resize(offset + ZSTD_compressBound(srcBytes));
  size_t compressed = ZSTD_compress(&out[offset], out.size() - offset, src, srcBytes, kZstdLevel);
  if (ZSTD_isError(compressed)) {
    // Writing no data makes readers reject this record instead of misreading it.
    SpzLog("[SPZ ERROR] Custom attribute '%s': zstd compression failed", name.c_str());
    compressed = 0;
  }
  out.resize(offset + compressed);
  return out;
}

uint32_t SpzExtensionCustomAttributeAdobe::payloadBytes() const {
  return static_cast<uint32_t>(encode().size());
}

void SpzExtensionCustomAttributeAdobe::write(std::ostream& os) const {
  const std::string payload = encode();
  const uint32_t t = static_cast<uint32_t>(extensionType);
  const uint32_t len = static_cast<uint32_t>(payload.size());  // validate() keeps this under 4 GiB
  SpzLog("[SPZ] Writing extension: CustomAttributeAdobe (name=%s, dim=%u, bits=%u, %u bytes)",
         name.c_str(), static_cast<unsigned>(dim), static_cast<unsigned>(bits),
         static_cast<unsigned>(len));
  os.write(reinterpret_cast<const char*>(&t), sizeof(t));
  os.write(reinterpret_cast<const char*>(&len), sizeof(len));
  os.write(payload.data(), static_cast<std::streamsize>(payload.size()));
}

std::optional<SpzExtensionBasePtr> SpzExtensionCustomAttributeAdobe::read(
    std::istream& is, int32_t numPoints) {
  auto rec = std::make_shared<SpzExtensionCustomAttributeAdobe>();
  uint16_t nameLength{};
  uint8_t kindValue{};
  if (!readExact(is, nameLength)) {
    SpzLog("[SPZ ERROR] CustomAttributeAdobe: truncated payload");
    return std::nullopt;
  }
  rec->name.resize(nameLength);
  if (!is.read(rec->name.data(), nameLength) || !readExact(is, rec->dim) ||
      !readExact(is, kindValue) || !readExact(is, rec->bits)) {
    SpzLog("[SPZ ERROR] CustomAttributeAdobe: truncated payload");
    return std::nullopt;
  }
  rec->kind = static_cast<SpzCustomAttributeKind>(kindValue);
  const bool quantized = isQuantized(rec->bits);
  if (nameLength == 0 || rec->dim == 0 ||
      kindValue > static_cast<uint8_t>(SpzCustomAttributeKind::Vector3) ||
      (rec->kind == SpzCustomAttributeKind::Vector3 && rec->dim % 3 != 0) ||
      (!quantized && rec->bits != kFloatBits)) {
    SpzLog("[SPZ ERROR] CustomAttributeAdobe '%s': invalid settings (dim=%u, kind=%u, bits=%u)",
           rec->name.c_str(), static_cast<unsigned>(rec->dim), static_cast<unsigned>(kindValue),
           static_cast<unsigned>(rec->bits));
    return std::nullopt;
  }

  const size_t dim = rec->dim;
  std::vector<float> lo, hi;
  if (quantized) {
    lo.resize(dim);
    hi.resize(dim);
    if (!is.read(reinterpret_cast<char*>(lo.data()), dim * sizeof(float)) ||
        !is.read(reinterpret_cast<char*>(hi.data()), dim * sizeof(float))) {
      SpzLog("[SPZ ERROR] CustomAttributeAdobe '%s': truncated payload", rec->name.c_str());
      return std::nullopt;
    }
    for (size_t c = 0; c < dim; ++c) {
      if (!std::isfinite(lo[c]) || !std::isfinite(hi[c]) || lo[c] > hi[c]) {
        SpzLog("[SPZ ERROR] CustomAttributeAdobe '%s': invalid range", rec->name.c_str());
        return std::nullopt;
      }
    }
  }

  // The payload stream is an in-memory buffer, so in_avail() is the rest of it: the zstd frame.
  const std::streamsize compressedSize = is.rdbuf()->in_avail();
  const uint64_t expected =
      static_cast<uint64_t>(numPoints) * dim * (quantized ? 1 : sizeof(float));
  if (compressedSize <= 0 || expected > std::numeric_limits<size_t>::max() ||
      expected > static_cast<uint64_t>(compressedSize) * kMaxZstdExpansion) {
    SpzLog("[SPZ ERROR] CustomAttributeAdobe '%s': data doesn't cover %d points",
           rec->name.c_str(), numPoints);
    return std::nullopt;
  }
  std::string compressed(static_cast<size_t>(compressedSize), '\0');
  is.read(compressed.data(), compressedSize);

  const size_t points = static_cast<size_t>(numPoints);
  rec->values.resize(points * dim);
  std::vector<uint8_t> levels(quantized ? rec->values.size() : 0);
  std::vector<float> planar(!quantized && dim > 1 ? rec->values.size() : 0);
  void* dst = quantized ? static_cast<void*>(levels.data())
              : dim > 1 ? static_cast<void*>(planar.data())
                        : static_cast<void*>(rec->values.data());
  const size_t decompressed =
      ZSTD_decompress(dst, static_cast<size_t>(expected), compressed.data(), compressed.size());
  if (ZSTD_isError(decompressed) || decompressed != expected) {
    SpzLog("[SPZ ERROR] CustomAttributeAdobe '%s': data doesn't cover %d points",
           rec->name.c_str(), numPoints);
    return std::nullopt;
  }

  if (quantized) {
    const uint8_t maxLevel = static_cast<uint8_t>((1 << rec->bits) - 1);
    for (size_t c = 0; c < dim; ++c) {
      const double step = (static_cast<double>(hi[c]) - lo[c]) / maxLevel;
      for (size_t p = 0; p < points; ++p) {
        const uint8_t level = levels[c * points + p];
        if (level > maxLevel) {
          SpzLog("[SPZ ERROR] CustomAttributeAdobe '%s': level out of range", rec->name.c_str());
          return std::nullopt;
        }
        rec->values[p * dim + c] = static_cast<float>(lo[c] + level * step);
      }
    }
  } else if (dim > 1) {
    for (size_t c = 0; c < dim; ++c) {
      for (size_t p = 0; p < points; ++p)
        rec->values[p * dim + c] = planar[c * points + p];
    }
  }

  SpzLog("[SPZ] Found extension: CustomAttributeAdobe (name=%s, dim=%u, bits=%u)",
         rec->name.c_str(), static_cast<unsigned>(rec->dim), static_cast<unsigned>(rec->bits));
  return std::optional{std::move(rec)};
}

SpzExtensionType SpzExtensionCustomAttributeAdobe::type() {
  return SpzExtensionType::SPZ_ADOBE_custom_attribute;
}

SpzExtensionBase* SpzExtensionCustomAttributeAdobe::copyAsRawData() const {
  return new SpzExtensionCustomAttributeAdobe(*this);
}

// No PLY representation: custom attributes would be vertex properties, which PLY extensions can't add.
std::optional<SpzExtensionBasePtr> SpzExtensionCustomAttributeAdobe::tryReadFromPly(
    std::istream&, const std::vector<PlyExtraElement>&) const {
  return std::nullopt;
}

void SpzExtensionCustomAttributeAdobe::writePlyHeader(std::ostream&) const {}

void SpzExtensionCustomAttributeAdobe::writePlyData(std::ostream&) const {}

}  // namespace spz
