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

}  // namespace spz
