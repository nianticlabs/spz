/*
MIT License

Copyright (c) 2025 Niantic Labs
Copyright (c) 2025 Adobe Inc.

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

#pragma once
#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "splat-types.h"
#include "splat-extensions.h"

namespace spz {
#ifdef ANDROID
static constexpr char LOG_TAG[] = "SPZ";
template <class... Args>
static void SpzLog(const char *fmt, Args &&...args) {
  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, fmt, std::forward<Args>(args)...);
}
#else
template <class... Args>
static void SpzLog(const char *fmt, Args &&...args) {
  printf(fmt, std::forward<Args>(args)...);
  printf("\n");
  fflush(stdout);
}
#endif  // ANDROID

template <class... Args>
static void SpzLog(const char *fmt) {
  SpzLog("%s", fmt);
}

// Maximum degree supported
constexpr int SH_MAX_DEGREE = 4;
constexpr int SH_MAX_COEFFS = 24;
constexpr int DEFAULT_SH1_BITS = 5;
constexpr int DEFAULT_SH_REST_BITS = 4;

// Latest version of the packed format, update this when changing the format.
constexpr int LATEST_SPZ_HEADER_VERSION = 2;

// Represents a single inflated gaussian. Each gaussian has 236 bytes. Although the data is easier
// to intepret in this format, it is not more precise than the packed format, since it was inflated.
struct UnpackedGaussian {
  std::array<float, 3> position;  // x, y, z
  std::array<float, 4> rotation;  // x, y, z, w
  std::array<float, 3> scale;     // std::log(scale)
  std::array<float, 3> color;     // rgb sh0 encoding
  float alpha;                    // inverse logistic
  std::array<float, SH_MAX_COEFFS> shR;
  std::array<float, SH_MAX_COEFFS> shG;
  std::array<float, SH_MAX_COEFFS> shB;
};

// Represents a single low precision gaussian. Each gaussian has exactly 64 bytes, even if it does
// not have full spherical harmonics.
struct PackedGaussian {
  std::array<uint8_t, 9> position{};
  std::array<uint8_t, 3> rotation{};
  std::array<uint8_t, 3> scale{};
  std::array<uint8_t, 3> color{};
  uint8_t alpha = 0;
  std::array<uint8_t, SH_MAX_COEFFS> shR{};
  std::array<uint8_t, SH_MAX_COEFFS> shG{};
  std::array<uint8_t, SH_MAX_COEFFS> shB{};

  UnpackedGaussian unpack(
    bool usesFloat16, int32_t fractionalBits, const CoordinateConverter &c,
    float shMin = -1.0f, float shMax = 1.0f) const;
};

// Represents a full splat with lower precision. Each splat has at most 64 bytes, although splats
// with fewer spherical harmonics degrees will have less. The data is stored non-interleaved.
struct PackedGaussians {
  uint32_t version = LATEST_SPZ_HEADER_VERSION;  // Version of the packed format
  int32_t numPoints = 0;       // Total number of points (gaussians)
  int32_t shDegree = 0;        // Degree of spherical harmonics
  int32_t fractionalBits = 0;  // Number of bits used for fractional part of fixed-point coords
  bool antialiased = false;    // Whether gaussians should be rendered with mip-splat antialiasing

  std::vector<uint8_t> positions;
  std::vector<uint8_t> scales;
  std::vector<uint8_t> rotations;
  std::vector<uint8_t> alphas;
  std::vector<uint8_t> colors;
  std::vector<uint8_t> sh;

  std::vector<SpzExtensionBasePtr> extensions;  // List of extensions, if any

  bool usesFloat16() const;
  PackedGaussian at(int32_t i) const;
  UnpackedGaussian unpack(int32_t i, const CoordinateConverter &c) const;
};

struct PackOptions {
  uint32_t version = LATEST_SPZ_HEADER_VERSION;  // Version of the packed format

  CoordinateSystem from = CoordinateSystem::UNSPECIFIED;
  // Spherical harmonics quantization parameters
  uint8_t sh1Bits = DEFAULT_SH1_BITS;      // Bits for SH degree 1 coefficients (max 8)
  uint8_t shRestBits = DEFAULT_SH_REST_BITS;   // Bits for SH degree 2+ coefficients (max 8)
  bool enableSHMinMaxScaling = false;  // Whether to enable SH min/max scaling, useful for legacy versions.
};

struct UnpackOptions {
  CoordinateSystem to = CoordinateSystem::UNSPECIFIED;
};

// Saves Gaussian splat in packed format, returning a vector of bytes.
bool saveSpz(
  const GaussianCloud &gaussians, const PackOptions &options, std::vector<uint8_t> *output);

// Loads Gaussian splat from a vector of bytes in packed format.
GaussianCloud loadSpz(const std::vector<uint8_t> &data, const UnpackOptions &options);

// Loads Gaussian splat from a vector of bytes in packed format.
PackedGaussians loadSpzPacked(const std::string &filename);
PackedGaussians loadSpzPacked(const uint8_t *data, int32_t size);
PackedGaussians loadSpzPacked(const std::vector<uint8_t> &data);

// Saves Gaussian splat in packed format to a file
bool saveSpz(
  const GaussianCloud &gaussians, const PackOptions &options, const std::string &filename);

// Loads Gaussian splat from a file in packed format
GaussianCloud loadSpz(const std::string &filename, const UnpackOptions &o);

// Saves Gaussian splat data in .ply format
bool saveSplatToPly(
  const spz::GaussianCloud &gaussians, const PackOptions &options, const std::string &filename);

// Loads Gaussian splat data in .ply format
GaussianCloud loadSplatFromPly(const std::string &filename, const UnpackOptions &options);

void serializePackedGaussians(const PackedGaussians &packed, std::ostream *out);

bool compressGzipped(const uint8_t *data, size_t size, std::vector<uint8_t> *out);
}  // namespace spz
