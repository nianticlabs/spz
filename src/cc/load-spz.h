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
#ifdef SPZ_BUILD_EXTENSIONS
#include "splat-extensions.h"
#endif

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

// Number of bits used for the fractional part of fixed-point position
// coordinates. Positions are stored as 24-bit signed fixed-point per axis, so
// the representable range is +/- 2^(23 - fractionalBits). With the default 12,
// this is ~0.25 mm resolution over a +/- 2048 m range.
constexpr int DEFAULT_FRACTIONAL_BITS = 12;

// Maximum number of fractional bits that fit in the 24-bit signed integer
// representation of positions. At this maximum, all 23 value bits are
// fractional and the per-axis representable range collapses to +/- 1.0.
// Higher values would either silently truncate the encoded fixed-point value
// (it would no longer fit in int24) or invoke undefined behavior in the
// (1 << fractionalBits) computation used by both encode and decode.
constexpr int MAX_FRACTIONAL_BITS = 23;

// Minimum number of fractional bits accepted by the encoder. Below this the
// resolution is too coarse for any realistic Gaussian-splat use case: at N=4
// the per-axis range is +/- 524288 m and resolution is 6.25 cm, already
// generous for aerial / city-scale scenes. N=3 has 12.5 cm resolution over a
// 1048 km range, and the lower values get progressively worse -- almost
// certainly a configuration mistake, so reject loudly instead.
constexpr int MIN_FRACTIONAL_BITS = 4;

// Latest version of the packed format, update this when changing the format.
constexpr int LATEST_SPZ_HEADER_VERSION = 4;

// Minimum version of the ZSTD-compressed SPZ format.
constexpr int MIN_ZSTD_SPZ_HEADER_VERSION = 4;

// Minimum version of that uses SmallestThree quaternions.
constexpr int MIN_SMALLEST_THREE_QUATERNIONS_VERSION = 3;

// NGSP Magic Number used for spz file identification.
constexpr uint32_t NGSP_MAGIC = 0x5053474e;

// Represents a single inflated gaussian. Each gaussian has 344 bytes (position, rotation, scale,
// color, alpha, and 24 SH coeffs x 3 channels). Although the data is easier to interpret in this
// format, it is not more precise than the packed format, since it was inflated.
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

// Represents a single low precision gaussian. Each gaussian has exactly 92 bytes (for degree 4
// spherical harmonics); the struct layout is fixed so that at() can index into non-interleaved buffers.
struct PackedGaussian {
  std::array<uint8_t, 9> position{};
  std::array<uint8_t, 4> rotation{};
  std::array<uint8_t, 3> scale{};
  std::array<uint8_t, 3> color{};
  uint8_t alpha = 0;
  std::array<uint8_t, SH_MAX_COEFFS> shR{};
  std::array<uint8_t, SH_MAX_COEFFS> shG{};
  std::array<uint8_t, SH_MAX_COEFFS> shB{};

  UnpackedGaussian unpack(
    bool usesFloat16, bool usesQuaternionSmallestThree, int32_t fractionalBits, const CoordinateConverter &c) const;
};

// Represents a full splat with lower precision. Each splat has at most 92 bytes (for degree 4
// spherical harmonics), although splats with fewer spherical harmonics degrees will have less.
// The data is stored non-interleaved.
struct PackedGaussians {
  uint32_t version = LATEST_SPZ_HEADER_VERSION;  // Version of the packed format
  int32_t numPoints = 0;       // Total number of points (gaussians)
  int32_t shDegree = 0;        // Degree of spherical harmonics
  int32_t fractionalBits = 0;  // Number of bits used for fractional part of fixed-point coords
  bool antialiased = false;    // Whether gaussians should be rendered with mip-splat antialiasing
  bool usesQuaternionSmallestThree = true; // Whether gaussians use the smallest three method to store quaternions

  std::vector<uint8_t> positions;
  std::vector<uint8_t> scales;
  std::vector<uint8_t> rotations;
  std::vector<uint8_t> alphas;
  std::vector<uint8_t> colors;
  std::vector<uint8_t> sh;

#ifdef SPZ_BUILD_EXTENSIONS
  std::vector<SpzExtensionBasePtr> extensions;  // List of extensions, if any
#endif

  bool usesFloat16() const;
  PackedGaussian at(int32_t i) const;
  UnpackedGaussian unpack(int32_t i, const CoordinateConverter &c) const;
};

struct PackOptions {
  uint32_t version = LATEST_SPZ_HEADER_VERSION;  // Version of the packed format

  CoordinateSystem from = CoordinateSystem::UNSPECIFIED;

  // Quantization bits are only used during packing to reduce information entropy for g-zipping.
  // Unpacking doesn't need these values since g-unzipping already fills zero bits for quantized data.
  uint8_t sh1Bits = DEFAULT_SH1_BITS;     // Bits for SH degree 1 coefficients
  uint8_t shRestBits = DEFAULT_SH_REST_BITS;  // Bits for SH degree 2+ coefficients

  // Number of fractional bits for fixed-point position coordinates. See
  // DEFAULT_FRACTIONAL_BITS for the meaning. Must satisfy
  // MIN_FRACTIONAL_BITS <= fractionalBits <= MAX_FRACTIONAL_BITS. The decoder
  // reads this value from the file header, so any value in range is fully
  // backward-compatible.
  uint8_t fractionalBits = DEFAULT_FRACTIONAL_BITS;
};

struct UnpackOptions {
  CoordinateSystem to = CoordinateSystem::UNSPECIFIED;
};

// Structure for PLY extra elements (non-vertex elements)
struct PlyExtraElement {
  std::string name;
  int32_t count;
  size_t bytesPerElement;
  bool isKnown;  // true for elements we explicitly handle (like safe_orbit)
};

// Saves Gaussian splat in packed format, returning a vector of bytes.
bool saveSpz(
  const GaussianCloud &gaussians, const PackOptions &options, std::vector<uint8_t> *output);

// Loads Gaussian splat from a vector of bytes in packed format.
GaussianCloud loadSpz(const std::vector<uint8_t> &data, const UnpackOptions &options);

// Loads Gaussian splat from a file / byte pointer / vector in packed format.
PackedGaussians loadSpzPacked(const std::string &filename);
PackedGaussians loadSpzPacked(const uint8_t *data, size_t size);
PackedGaussians loadSpzPacked(const std::vector<uint8_t> &data);

// Saves Gaussian splat in packed format to a file
bool saveSpz(
  const GaussianCloud &gaussians, const PackOptions &options, const std::string &filename);

// Loads Gaussian splat from a file in packed format
GaussianCloud loadSpz(const std::string &filename, const UnpackOptions &o);

// Loads Gaussian splat from a byte pointer in packed format.
GaussianCloud loadSpz(const uint8_t *data, size_t size, const UnpackOptions &options);

// Saves Gaussian splat data in .ply format
bool saveSplatToPly(
  const spz::GaussianCloud &gaussians, const PackOptions &options, const std::string &filename);

// Loads Gaussian splat data in .ply format
GaussianCloud loadSplatFromPly(const std::string &filename, const UnpackOptions &options);

void serializePackedGaussians(const PackedGaussians &packed, std::ostream *out);

bool compressGzipped(const uint8_t *data, size_t size, std::vector<uint8_t> *out);

// Returns true if the build has extension support enabled, false otherwise
bool hasExtensionSupport();
}  // namespace spz
