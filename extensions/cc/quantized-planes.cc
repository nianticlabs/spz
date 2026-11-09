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

#include "quantized-planes.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace spz {
namespace {
constexpr uint8_t kBitsPerByte = 8;

uint32_t maxLevel(uint8_t bits) { return (1u << bits) - 1; }

// Planes is a template parameter so the per-byte loops unroll; with one plane they vanish.
template <size_t Planes>
void encode(const float* values, size_t stride, size_t points, double min, double scale,
            double top, uint8_t* out) {
  for (size_t p = 0; p < points; ++p) {
    const double level = std::clamp(std::round((values[p * stride] - min) * scale), 0.0, top);
    const auto levelBits = static_cast<uint32_t>(level);
    for (size_t k = 0; k < Planes; ++k) {
      out[k * points + p] = static_cast<uint8_t>(levelBits >> (k * kBitsPerByte));
    }
  }
}

template <size_t Planes>
bool decode(const uint8_t* in, size_t points, double min, double step, uint32_t top,
            float* values, size_t stride) {
  for (size_t p = 0; p < points; ++p) {
    uint32_t level = 0;
    for (size_t k = 0; k < Planes; ++k) {
      level |= static_cast<uint32_t>(in[k * points + p]) << (k * kBitsPerByte);
    }
    if (level > top) return false;
    values[p * stride] = static_cast<float>(min + level * step);
  }
  return true;
}
}  // namespace

size_t quantizedLevelBytes(uint8_t bits) { return (bits + kBitsPerByte - 1) / kBitsPerByte; }

void encodeQuantizedPlanes(const float* values, size_t stride, size_t points, float min, float max,
                           uint8_t bits, uint8_t* out) {
  const double top = maxLevel(bits);
  // Double precision keeps the range finite even when max - min overflows float.
  const double range = static_cast<double>(max) - min;
  const double scale = range > 0.0 ? top / range : 0.0;
  if (quantizedLevelBytes(bits) == 1) {
    encode<1>(values, stride, points, min, scale, top, out);
  } else {
    encode<2>(values, stride, points, min, scale, top, out);
  }
}

bool decodeQuantizedPlanes(const uint8_t* in, size_t points, float min, float max, uint8_t bits,
                           float* values, size_t stride) {
  const uint32_t top = maxLevel(bits);
  const double step = (static_cast<double>(max) - min) / top;
  return quantizedLevelBytes(bits) == 1 ? decode<1>(in, points, min, step, top, values, stride)
                                        : decode<2>(in, points, min, step, top, values, stride);
}

std::optional<size_t> componentMajorBytes(size_t points, uint16_t dim, uint8_t bits) {
  const size_t perComponent = isQuantizedBits(bits) ? quantizedLevelBytes(bits) : sizeof(float);
  // points < 2^32, dim < 2^16 and perComponent <= 4, so this fits in 64 bits.
  const uint64_t bytes = static_cast<uint64_t>(points) * dim * perComponent;
  if (bytes > std::numeric_limits<size_t>::max()) return std::nullopt;
  return static_cast<size_t>(bytes);
}

void componentRanges(const std::vector<float>& values, uint16_t dim, std::vector<float>* min,
                     std::vector<float>* max) {
  min->assign(dim, 0.0f);
  max->assign(dim, 0.0f);
  const size_t points = dim == 0 ? 0 : values.size() / dim;
  if (points == 0) return;
  std::copy(values.begin(), values.begin() + dim, min->begin());
  std::copy(values.begin(), values.begin() + dim, max->begin());
  for (size_t p = 1; p < points; ++p) {
    for (size_t c = 0; c < dim; ++c) {
      (*min)[c] = std::min((*min)[c], values[p * dim + c]);
      (*max)[c] = std::max((*max)[c], values[p * dim + c]);
    }
  }
}

// Byte planes group similar bytes together, which compresses better than whole levels.
std::vector<uint8_t> encodeComponentMajor(const std::vector<float>& values, size_t points,
                                          uint16_t dim, uint8_t bits, const std::vector<float>& min,
                                          const std::vector<float>& max) {
  std::vector<uint8_t> out(*componentMajorBytes(points, dim, bits));
  if (points == 0) return out;
  if (isQuantizedBits(bits)) {
    const size_t componentBytes = points * quantizedLevelBytes(bits);
    for (size_t c = 0; c < dim; ++c) {
      encodeQuantizedPlanes(&values[c], dim, points, min[c], max[c], bits, &out[c * componentBytes]);
    }
  } else {
    for (size_t c = 0; c < dim; ++c) {
      for (size_t p = 0; p < points; ++p) {
        std::memcpy(&out[(c * points + p) * sizeof(float)], &values[p * dim + c], sizeof(float));
      }
    }
  }
  return out;
}

bool decodeComponentMajor(const uint8_t* bytes, size_t points, uint16_t dim, uint8_t bits,
                          const std::vector<float>& min, const std::vector<float>& max,
                          std::vector<float>* values) {
  values->resize(points * dim);
  if (points == 0) return true;
  if (isQuantizedBits(bits)) {
    const size_t componentBytes = points * quantizedLevelBytes(bits);
    for (size_t c = 0; c < dim; ++c) {
      if (!decodeQuantizedPlanes(bytes + c * componentBytes, points, min[c], max[c], bits,
                                 &(*values)[c], dim)) {
        return false;
      }
    }
  } else {
    for (size_t c = 0; c < dim; ++c) {
      for (size_t p = 0; p < points; ++p) {
        std::memcpy(&(*values)[p * dim + c], bytes + (c * points + p) * sizeof(float),
                    sizeof(float));
      }
    }
  }
  return true;
}

namespace {
// Signed levels zigzag to unsigned ones so small magnitudes of either sign keep the high bytes 0.
uint32_t zigzag(int32_t level) {
  return (static_cast<uint32_t>(level) << 1) ^ static_cast<uint32_t>(level >> 31);
}

int32_t unzigzag(uint32_t z) {
  return static_cast<int32_t>(z >> 1) ^ -static_cast<int32_t>(z & 1);
}

// The level for v, or std::nullopt if v isn't finite or the level doesn't fit in 32 bits.
std::optional<int32_t> stepLevel(float v, float step) {
  const double level = std::round(static_cast<double>(v) / step);
  if (!std::isfinite(level) || level < std::numeric_limits<int32_t>::min() ||
      level > std::numeric_limits<int32_t>::max()) {
    return std::nullopt;
  }
  return static_cast<int32_t>(level);
}

float levelValue(int32_t level, float step) {
  return static_cast<float>(static_cast<double>(level) * step);
}
}  // namespace

bool validSteps(const std::vector<float>& steps) {
  return std::all_of(steps.begin(), steps.end(),
                     [](float s) { return std::isfinite(s) && s > 0.0f; });
}

bool encodeComponentMajorSteps(const std::vector<float>& values, size_t points, uint16_t dim,
                               const std::vector<float>& steps, std::vector<uint8_t>* out) {
  out->assign(points * dim * kStepLevelBytes, 0);
  for (size_t c = 0; c < dim; ++c) {
    uint8_t* planes = out->data() + c * points * kStepLevelBytes;
    for (size_t p = 0; p < points; ++p) {
      const std::optional<int32_t> level = stepLevel(values[p * dim + c], steps[c]);
      if (!level) return false;
      const uint32_t z = zigzag(*level);
      for (size_t k = 0; k < kStepLevelBytes; ++k) {
        planes[k * points + p] = static_cast<uint8_t>(z >> (k * kBitsPerByte));
      }
    }
  }
  return true;
}

void decodeComponentMajorSteps(const uint8_t* bytes, size_t points, uint16_t dim,
                               const std::vector<float>& steps, std::vector<float>* values) {
  values->resize(points * dim);
  for (size_t c = 0; c < dim; ++c) {
    const uint8_t* planes = bytes + c * points * kStepLevelBytes;
    for (size_t p = 0; p < points; ++p) {
      uint32_t z = 0;
      for (size_t k = 0; k < kStepLevelBytes; ++k) {
        z |= static_cast<uint32_t>(planes[k * points + p]) << (k * kBitsPerByte);
      }
      (*values)[p * dim + c] = levelValue(unzigzag(z), steps[c]);
    }
  }
}

float quantizeToStep(float v, float step) { return levelValue(*stepLevel(v, step), step); }

}  // namespace spz

