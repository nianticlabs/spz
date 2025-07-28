#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "splat-types.h"

namespace spz {

// Maximum degree supported
constexpr int SH_MAX_DEGREE = 4;
constexpr int SH_MAX_COEFFS = 24;

// Latest version of the packed format, update this when changing the format.
constexpr int LATEST_SPZ_HEADER_VERSION = 4;

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

  // SH quantization parameters
  uint8_t sh1Bits = 5;      // Bits for SH degree 1 coefficients
  uint8_t shRestBits = 4;   // Bits for SH degree 2+ coefficients
  float shMin = -1.0f;      // Minimum SH coefficient value used for quantization
  float shMax = 1.0f;       // Maximum SH coefficient value used for quantization

  // Safe orbit camera parameters
  bool hasSafeOrbit = false;           // Whether safe orbit data is present
  float safeOrbitElevationMin = 0.0f;  // Minimum elevation for safe orbit (radians)
  float safeOrbitElevationMax = 0.0f;  // Maximum elevation for safe orbit (radians)
  float safeOrbitRadiusMin = 0.0f;     // Minimum radius for safe orbit

  std::vector<uint8_t> positions;
  std::vector<uint8_t> scales;
  std::vector<uint8_t> rotations;
  std::vector<uint8_t> alphas;
  std::vector<uint8_t> colors;
  std::vector<uint8_t> sh;

  bool usesFloat16() const;
  PackedGaussian at(int32_t i) const;
  UnpackedGaussian unpack(int32_t i, const CoordinateConverter &c) const;
};

struct PackOptions {
  uint32_t version = LATEST_SPZ_HEADER_VERSION;  // Version of the packed format

  CoordinateSystem from = CoordinateSystem::UNSPECIFIED;
  // Spherical harmonics quantization parameters
  uint8_t sh1Bits = 5;      // Bits for SH degree 1 coefficients (max 8)
  uint8_t shRestBits = 4;   // Bits for SH degree 2+ coefficients (max 8)

  // Safe orbit camera parameters
  bool hasSafeOrbit = false;           // Whether safe orbit data is present
  float safeOrbitElevationMin = 0.0f;  // Minimum elevation for safe orbit (radians)
  float safeOrbitElevationMax = 0.0f;  // Maximum elevation for safe orbit (radians)
  float safeOrbitRadiusMin = 0.0f;     // Minimum radius for safe orbit
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
