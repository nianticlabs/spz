/*
MIT License

Copyright (c) 2026 Niantic Labs
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

#include "chunked-unpack.h"

#include <array>
#include <vector>

#include "splat-types.h"
#include "splat-utils.h"

namespace spz {

namespace {

void unpackPositions(const std::vector<uint8_t> &packed,
                     bool usesFloat16, int32_t fractionalBits,
                     const CoordinateConverter &c,
                     int32_t pointOffset, int32_t pointCount, float *out) {
  const int32_t coordCount = pointCount * 3;
  if (usesFloat16) {
    const auto *src =
        reinterpret_cast<const Half *>(packed.data()) + static_cast<size_t>(pointOffset) * 3;
    for (int32_t i = 0; i < coordCount; i++) {
      out[i] = halfToFloat(src[i]) * c.flipP[i % 3];
    }
  } else {
    const float scale = 1.0f / (1 << fractionalBits);
    const uint8_t *src = packed.data() + static_cast<size_t>(pointOffset) * 3 * 3;
    for (int32_t i = 0; i < coordCount; i++) {
      int32_t fixed32 = src[i * 3 + 0];
      fixed32 |= src[i * 3 + 1] << 8;
      fixed32 |= src[i * 3 + 2] << 16;
      fixed32 |= (fixed32 & 0x800000) ? 0xff000000 : 0;  // sign extension
      out[i] = static_cast<float>(fixed32) * scale * c.flipP[i % 3];
    }
  }
  if (c.rotFlipPFunc) {
    for (int32_t i = 0; i < pointCount; i++) {
      c.rotFlipPFunc(out + i * 3);
    }
  }
}

void unpackScales(const std::vector<uint8_t> &packed,
                  int32_t pointOffset, int32_t pointCount, float *out) {
  const int32_t coordCount = pointCount * 3;
  const uint8_t *src = packed.data() + static_cast<size_t>(pointOffset) * 3;
  for (int32_t i = 0; i < coordCount; i++) {
    out[i] = src[i] / 16.0f - 10.0f;
  }
}

void unpackRotations(const std::vector<uint8_t> &packed,
                     bool usesQuaternionSmallestThree,
                     const CoordinateConverter &c,
                     int32_t pointOffset, int32_t pointCount, float *out) {
  if (usesQuaternionSmallestThree) {
    const uint8_t *src = packed.data() + static_cast<size_t>(pointOffset) * 4;
    for (int32_t i = 0; i < pointCount; i++) {
      unpackQuaternionSmallestThree(&out[i * 4], &src[i * 4], c);
    }
  } else {
    const uint8_t *src = packed.data() + static_cast<size_t>(pointOffset) * 3;
    for (int32_t i = 0; i < pointCount; i++) {
      unpackQuaternionFirstThree(&out[i * 4], &src[i * 3], c);
    }
  }
}

void unpackAlphas(const std::vector<uint8_t> &packed,
                  int32_t pointOffset, int32_t pointCount, float *out) {
  const uint8_t *src = packed.data() + static_cast<size_t>(pointOffset);
  for (int32_t i = 0; i < pointCount; i++) {
    out[i] = invSigmoid(src[i] / 255.0f);
  }
}

void unpackColors(const std::vector<uint8_t> &packed,
                  int32_t pointOffset, int32_t pointCount, float *out) {
  const int32_t coordCount = pointCount * 3;
  const uint8_t *src = packed.data() + static_cast<size_t>(pointOffset) * 3;
  for (int32_t i = 0; i < coordCount; i++) {
    out[i] = ((src[i] / 255.0f) - 0.5f) / colorScale;
  }
}

void unpackSh(const std::vector<uint8_t> &packed, int32_t shDim,
              const CoordinateConverter &c,
              int32_t pointOffset, int32_t pointCount, float *out) {
  // Layout per point: shDim coefficients × 3 RGB channels. Within-family,
  // c.flipSh is the per-coefficient sign flip fused into unquantize below.
  const size_t fpp = static_cast<size_t>(shDim) * 3;
  const uint8_t *src = packed.data() + static_cast<size_t>(pointOffset) * fpp;
  for (int32_t p = 0; p < pointCount; p++) {
    const size_t base = static_cast<size_t>(p) * fpp;
    for (int32_t j = 0; j < shDim; j++) {
      const float flip = c.flipSh[j];
      out[base + j * 3 + 0] = flip * unquantizeSH(src[base + j * 3 + 0]);
      out[base + j * 3 + 1] = flip * unquantizeSH(src[base + j * 3 + 1]);
      out[base + j * 3 + 2] = flip * unquantizeSH(src[base + j * 3 + 2]);
    }
  }
  // Cross-family: coordinateConverter sets c.flipSh to all 1.0 (the loop above
  // is a no-op) and populates c.rotFlipShFuncs[band] with the analytic
  // rotation+flip per band. Mirrors the SH block in
  // GaussianCloud::convertCoordinates (splat-types.h) — if this grows a third
  // caller, factor the rotate-then-scatter into a shared helper.
  if (!c.rotFlipShFuncs[0]) {
    return;
  }
  for (int32_t p = 0; p < pointCount; p++) {
    const size_t base = static_cast<size_t>(p) * fpp;
    for (int band = 0; band < SH_MAX_DEGREE; ++band) {
      if (!c.rotFlipShFuncs[static_cast<size_t>(band)]) { break; }
      const size_t bandStart = static_cast<size_t>(band * (band + 2));
      const size_t bandSize = static_cast<size_t>(2 * band + 3);
      if (bandStart + bandSize > static_cast<size_t>(shDim)) { break; }
      std::array<float, 9> tmp{};
      for (int channel = 0; channel < 3; ++channel) {
        for (size_t k = 0; k < bandSize; ++k) {
          tmp[k] = out[base + (bandStart + k) * 3 + static_cast<size_t>(channel)];
        }
        c.rotFlipShFuncs[static_cast<size_t>(band)](tmp.data());
        for (size_t k = 0; k < bandSize; ++k) {
          out[base + (bandStart + k) * 3 + static_cast<size_t>(channel)] = tmp[k];
        }
      }
    }
  }
}

}  // namespace

bool unpackChunk(PackedGaussians &packed,
                 SplatAttribute attr,
                 int32_t pointOffset,
                 int32_t pointCount,
                 float *out,
                 CoordinateSystem coordFrom,
                 CoordinateSystem coordTo) {
  if (pointOffset < 0 || pointCount < 0) return false;
  if (pointOffset + pointCount > packed.numPoints) return false;
  if (pointCount == 0) return true;
  if (out == nullptr) return false;

  const std::vector<uint8_t> &buf = packedBuffer(packed, attr);
  if (buf.empty() && packed.numPoints > 0) {
    // Attribute never present or already released.
    return false;
  }

  const CoordinateConverter c = coordinateConverter(coordFrom, coordTo, packed.shDegree);
  const int32_t shDim = dimForDegree(packed.shDegree);

  switch (attr) {
    case SplatAttribute::Positions:
      unpackPositions(buf, packed.usesFloat16(), packed.fractionalBits,
                      c, pointOffset, pointCount, out);
      return true;
    case SplatAttribute::Scales:
      unpackScales(buf, pointOffset, pointCount, out);
      return true;
    case SplatAttribute::Rotations:
      unpackRotations(buf, packed.usesQuaternionSmallestThree,
                      c, pointOffset, pointCount, out);
      return true;
    case SplatAttribute::Alphas:
      unpackAlphas(buf, pointOffset, pointCount, out);
      return true;
    case SplatAttribute::Colors:
      unpackColors(buf, pointOffset, pointCount, out);
      return true;
    case SplatAttribute::Sh:
      if (shDim == 0) return true;  // no SH data
      unpackSh(buf, shDim, c, pointOffset, pointCount, out);
      return true;
  }
  return false;
}

}  // namespace spz
