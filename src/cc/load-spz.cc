#include "load-spz.h"

#include <zlib.h>

#ifdef ANDROID
#include <android/log.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
// Verbosity level 1 = global stats, 2 = per-Gaussian results, 3 = debug.
#define VERBOSE 1
#ifdef VERBOSE
#include <cassert>
#endif

namespace spz {

namespace {

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

// Scale factor for DC color components. To convert to RGB, we should multiply by 0.282, but it can
// be useful to represent base colors that are out of range if the higher spherical harmonics bands
// bring them back into range so we multiply by a smaller value.
constexpr float colorScale = 0.15f;

int32_t degreeForDim(int32_t dim) {
  if (dim < 3)
    return 0;
  if (dim < 8)
    return 1;
  if (dim < 15)
    return 2;
  return 3;
}

int32_t dimForDegree(int32_t degree) {
  switch (degree) {
    case 0:
      return 0;
    case 1:
      return 3;
    case 2:
      return 8;
    case 3:
      return 15;
    default:
      SpzLog("[SPZ: ERROR] Unsupported SH degree: %d\n", degree);
      return 0;
  }
}

uint8_t toUint8(float x) { return static_cast<uint8_t>(std::clamp(std::round(x), 0.0f, 255.0f)); }

// Quantizes to 8 bits, the round to nearest bucket center. 0 always maps to a bucket center.
uint8_t quantizeSH(float x, int32_t bucketSize) {
  int32_t q = static_cast<int>(std::round(x * 128.0f) + 128.0f);
  q = (q + bucketSize / 2) / bucketSize * bucketSize;
  return static_cast<uint8_t>(std::clamp(q, 0, 255));
}

float unquantizeSH(uint8_t x) { return (static_cast<float>(x) - 128.0f) / 128.0f; }

float sigmoid(float x) { return 1 / (1 + std::exp(-x)); }

float invSigmoid(float x) { return std::log(x / (1.0f - x)); }

template <typename T>
size_t countBytes(std::vector<T> vec) {
  return vec.size() * sizeof(vec[0]);
}

#define CHECK(x)                                                              \
  {                                                                           \
    if (!(x)) {                                                               \
      SpzLog("[SPZ: ERROR] Check failed: %s:%d: %s", __FILE__, __LINE__, #x); \
      return false;                                                           \
    }                                                                         \
  }

#define CHECK_GE(x, y) CHECK((x) >= (y))
#define CHECK_LE(x, y) CHECK((x) <= (y))
#define CHECK_EQ(x, y) CHECK((x) == (y))

bool checkSizes(const GaussianCloud &g) {
  CHECK_GE(g.numPoints, 0);
  CHECK_GE(g.shDegree, 0);
  CHECK_LE(g.shDegree, 3);
  CHECK_EQ(g.positions.size(), g.numPoints * 3);
  CHECK_EQ(g.scales.size(), g.numPoints * 3);
  CHECK_EQ(g.rotations.size(), g.numPoints * 4);
  CHECK_EQ(g.alphas.size(), g.numPoints);
  CHECK_EQ(g.colors.size(), g.numPoints * 3);
  CHECK_EQ(g.sh.size(), g.numPoints * dimForDegree(g.shDegree) * 3);
  return true;
}

bool checkSizes(const PackedGaussians &packed, int32_t numPoints, int32_t shDim, bool usesFloat16) {
  CHECK_EQ(packed.positions.size(), numPoints * 3 * (usesFloat16 ? 2 : 3));
  CHECK_EQ(packed.scales.size(), numPoints * 3);
  CHECK_EQ(packed.rotations.size(), numPoints * 3);
  CHECK_EQ(packed.alphas.size(), numPoints);
  CHECK_EQ(packed.colors.size(), numPoints * 3);
  CHECK_EQ(packed.sh.size(), numPoints * shDim * 3);
  return true;
}

constexpr uint8_t FlagAntialiased = 0x1;

struct PackedGaussiansHeader {
  uint32_t magic = 0x5053474e;  // NGSP = Niantic gaussian splat
  uint32_t version = 2;
  uint32_t numPoints = 0;
  uint8_t shDegree = 0;
  uint8_t fractionalBits = 0;
  uint8_t flags = 0;
  uint8_t reserved = 0;
};

bool decompressGzippedImpl(
  const uint8_t *compressed, size_t size, int32_t windowSize, std::vector<uint8_t> *out) {
  std::vector<uint8_t> buffer(8192);
  z_stream stream = {};
  stream.next_in = const_cast<Bytef *>(compressed);
  stream.avail_in = size;
  if (inflateInit2(&stream, windowSize) != Z_OK) {
    return false;
  }
  out->clear();
  bool success = false;
  while (true) {
    stream.next_out = buffer.data();
    stream.avail_out = buffer.size();
    int32_t res = inflate(&stream, Z_NO_FLUSH);
    if (res != Z_OK && res != Z_STREAM_END) {
      break;
    }
    out->insert(out->end(), buffer.data(), buffer.data() + buffer.size() - stream.avail_out);
    if (res == Z_STREAM_END) {
      success = true;
      break;
    }
  }
  inflateEnd(&stream);
  return success;
}

bool decompressGzipped(const uint8_t *compressed, size_t size, std::vector<uint8_t> *out) {
  // Here 16 means enable automatic gzip header detection; consider switching this to 32 to enable
  // both automated gzip and zlib header detection.
  return decompressGzippedImpl(compressed, size, 16 | MAX_WBITS, out);
}

bool decompressGzipped(const uint8_t *compressed, size_t size, std::string *out) {
  std::vector<uint8_t> buffer;
  if (!decompressGzipped(compressed, size, &buffer)) {
    return false;
  }
  out->assign(reinterpret_cast<const char *>(buffer.data()), buffer.size());
  return true;
}

// unpack an 8-bit quantized quaternion
// returns a normalized quaternion (unit norm with s[3]>=0)
Quat4f unpackQuaternion(const uint8_t* r) {
    // unpack
    const Vec3f xyz = plus(
        times(
            Vec3f{static_cast<float>(r[0]), static_cast<float>(r[1]), static_cast<float>(r[2])},
            1.0f / 127.5f),
        Vec3f{-1, -1, -1});
    Quat4f s;
    std::copy(xyz.data(), xyz.data() + 3, &s[0]);
    // Compute the real component - we know the quaternion is normalized and w is non-negative
    float xyz_squaredNorm = squaredNorm(xyz);
    if (xyz_squaredNorm <= 1) {
        s[3] = std::sqrt(1.0f - squaredNorm(xyz));
    } else {
        s[3] = 0.f;
        s = times(s, 1.f / std::sqrt(xyz_squaredNorm));
    }
    return s;
}

// Compute the quaternion reconstruction error: 1 - dot(q,r)
// q and r are supposed to be unit quaternions, and their scalar component
// (4th component) should be positive.
// Ref: https://math.stackexchange.com/a/90098
float quaternionReconstructionError(const Quat4f& q, const Quat4f& r) {
  return 1. - std::pow(dot(q, r), 2);
}

// Compute the geodesic distance in degrees: acos(clamp(2 * dot(q,r)^2 - 1.0, -1.0, 1.0))
// q and r are supposed to be unit quaternions, and their scalar component
// (4th component) should be positive.
// Note: 2 * acos(abs(dot(q,r))) would work too.
// Ref: https://math.stackexchange.com/a/90098
float quaternionGeodesicDistanceDegree(const Quat4f& q, const Quat4f& r) {
  // Compute the real component - we know the quaternion is normalized and w is non-negative
  return std::acos(std::clamp(2 * std::pow(dot(q, r), 2) - 1.0, -1.0, 1.0)) * 180 / M_PI;
}

constexpr float quaternionGeodesicDistanceDegreeToReconstructionError(float deg) {
  return (1 - std::cos(deg * M_PI / 180)) / 2;
}

// Compute the rotation angle distance in degrees: abs(2 * acos(q[3]) - 2.0 * std::acos(std::clamp(std::abs(s[3]), -1.0f, 1.0f)))
// q and r are supposed to be unit quaternions, and their scalar component
// (4th component) should be positive.
float quaternionRotationAngleDistanceDegree(const Quat4f& q, const Quat4f& r) {
  // For a unit quaternion [x,y,z,w], the rotation angle is 2*arccos(|w|)
  return std::abs(2.0 * std::acos(std::clamp(std::abs(q[3]), -1.0f, 1.0f)) - 2.0 * std::acos(std::clamp(std::abs(r[3]), -1.0f, 1.0f))) * 180 / M_PI;
}

// Compute the angle between rotation axes in degrees
// q and r are supposed to be unit quaternions, and their scalar component
// (4th component) should be positive.
float quaternionAxisAngleDistanceDegree(const Quat4f& q, const Quat4f& r) {
  // Get rotation axis direction for each quaternion
  Vec3f axis_q = normalized(Vec3f{q[0], q[1], q[2]});
  Vec3f axis_r = normalized(Vec3f{r[0], r[1], r[2]});

  // Compute angle between rotation axes in radians
  double dist_axis_angle_deg = std::acos(std::clamp(dot(axis_q, axis_r), -1.0f, 1.0f)) * 180 / M_PI;
  // Handle case where axes point in opposite directions but represent same rotation
  if (dist_axis_angle_deg > 90.0) {
    dist_axis_angle_deg = 180.0 - dist_axis_angle_deg;
  }
  return dist_axis_angle_deg;
}

}  // namespace

bool compressGzipped(const uint8_t *data, size_t size, std::vector<uint8_t> *out) {
  std::vector<uint8_t> buffer(8192);
  z_stream stream = {};
  if (
    deflateInit2(&stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 16 + MAX_WBITS, 9, Z_DEFAULT_STRATEGY)
    != Z_OK) {
    return false;
  }
  out->clear();
  out->reserve(size / 4);
  stream.next_in = const_cast<Bytef *>(reinterpret_cast<const Bytef *>(data));
  stream.avail_in = size;
  bool success = false;
  while (true) {
    stream.next_out = buffer.data();
    stream.avail_out = buffer.size();
    int32_t res = deflate(&stream, Z_FINISH);
    if (res != Z_OK && res != Z_STREAM_END) {
      break;
    }
    out->insert(out->end(), buffer.data(), buffer.data() + buffer.size() - stream.avail_out);
    if (res == Z_STREAM_END) {
      success = true;
      break;
    }
  }
  deflateEnd(&stream);
  return success;
}

PackedGaussians packGaussians(const GaussianCloud &g, const PackOptions &o) {
  if (!checkSizes(g)) {
    return {};
  }
  const int32_t numPoints = g.numPoints;
  const int32_t shDim = dimForDegree(g.shDegree);
  CoordinateConverter c = coordinateConverter(o.from, CoordinateSystem::RUB);

  // Use 12 bits for the fractional part of coordinates (~0.25 millimeter resolution). In the future
  // we can use different values on a per-splat basis and still be compatible with the decoder.
  PackedGaussians packed;
  packed.numPoints = g.numPoints;
  packed.shDegree = g.shDegree;
  packed.fractionalBits = 12;
  packed.antialiased = g.antialiased;
  packed.positions.resize(numPoints * 3 * 3);
  packed.scales.resize(numPoints * 3);
  packed.rotations.resize(numPoints * 3);
  packed.alphas.resize(numPoints);
  packed.colors.resize(numPoints * 3);
  packed.sh.resize(numPoints * shDim * 3);

  // Store coordinates as 24-bit fixed point values.
  const float scale = (1 << packed.fractionalBits);
  for (size_t i = 0; i < numPoints * 3; i++) {
    const int32_t fixed32 =
      static_cast<int32_t>(std::round(c.flipP[i % 3] * g.positions[i] * scale));
    packed.positions[i * 3 + 0] = fixed32 & 0xff;
    packed.positions[i * 3 + 1] = (fixed32 >> 8) & 0xff;
    packed.positions[i * 3 + 2] = (fixed32 >> 16) & 0xff;
  }

  for (size_t i = 0; i < numPoints * 3; i++) {
    packed.scales[i] = toUint8((g.scales[i] + 10.0f) * 16.0f);
  }

  // Define maximum search radius and error thresholds for early termination
  //
  // Maximum neighborhood size to explore (radius=3 means up to 7x7x7 cube)
  // Experimentally, it was determined that exploring beyond radius=8 never yields a better optimum.
  const int radius_max = 8;
  // Angle thresholds in degrees - if we find a representation with error below this, we can stop searching
  assert(radius_max >= 2);
  float geodesicDistanceThresholdDegree[radius_max] = {0.5, 1.0, 1.4};
  if(radius_max > 3) {
    geodesicDistanceThresholdDegree[3] = 1.7;
  }
  if(radius_max > 4) {
    geodesicDistanceThresholdDegree[4] = 2.1;
  }
  for (int radius = 5; radius < radius_max; ++radius) {
    const_cast<float*>(geodesicDistanceThresholdDegree)[radius] = 0.5f * radius;
  }

  // Convert angle thresholds to reconstruction error thresholds
  float errorThreshold[radius_max];
  for (int radius = 0; radius < radius_max; ++radius) {
    errorThreshold[radius] = quaternionGeodesicDistanceDegreeToReconstructionError(geodesicDistanceThresholdDegree[radius]);
#if VERBOSE >= 2
    if (!o.fast_rot_quantization) {
      SpzLog("[SPZ] error_threshold[%d]=%g", radius, errorThreshold[radius]);
    }
#endif
  }

  int32_t numPointsMaxRadius[radius_max + 1] = {0}; // Count of max radius used for testing points
  int32_t numPointsBestRadius[radius_max + 1] = {0}; // Count of best radius used for testing points
#if VERBOSE >= 1
  double averageGeodesicDistanceDegreeFast = 0.;
  double averageGeodesicDistanceDegreeBest = 0.;
#endif
  
  for (size_t i = 0; i < numPoints; i++) {
    // Normalize the quaternion, make w positive, then store xyz. w can be derived from xyz.
    // NOTE: These are already in xyzw order.
    Quat4f q = normalized(quat4f(&g.rotations[i * 4]));
  #if VERBOSE >= 2
    SpzLog("[SPZ] ** Gaussian %d", i );
    SpzLog("[SPZ] %g, %g, %g, %g (scalar)", g.rotations[i * 4 + 0], g.rotations[i * 4 + 1], g.rotations[i * 4 + 2], g.rotations[i * 4 + 3]);
    SpzLog("[SPZ] Normalized");
    SpzLog("[SPZ] %g, %g, %g, %g (scalar)", q[0], q[1], q[2], q[3]);
  #endif
    q[0] *= c.flipQ[0];
    q[1] *= c.flipQ[1];
    q[2] *= c.flipQ[2];

    Quat4f qq = times(q, (q[3] < 0 ? -127.5f : 127.5f));
    qq = plus(qq, Quat4f{127.5f, 127.5f, 127.5f, 127.5f});
    uint8_t r0[3] = {toUint8(qq[0]), toUint8(qq[1]), toUint8(qq[2])};

#if VERBOSE >= 1
    {
      // unpack
      Quat4f s = unpackQuaternion(r0);
      float dist = quaternionGeodesicDistanceDegree(q, s);
      averageGeodesicDistanceDegreeFast += dist;
#if VERBOSE >= 2
      SpzLog("[SPZ] Packed");
      SpzLog("[SPZ] %d, %d, %d", (int)r0[0], (int)r0[1], (int)r0[2]);
      // Distance between quaternions
      // Geodesic distance of the normalized quaternions. https://math.stackexchange.com/a/90098
      SpzLog("[SPZ] Unpacked normalized");
      SpzLog("[SPZ] %g, %g, %g, %g (scalar)", s[0], s[1], s[2], s[3]);
  
      SpzLog("[SPZ] Distances:");
      SpzLog("[SPZ] Reconstruction error: %g", quaternionGeodesicDistanceDegreeToReconstructionError(dist));
      SpzLog("[SPZ] Geodesic quaternion distance (degrees): %g", dist);
      SpzLog("[SPZ] Rotation angle distance (degrees): %g", quaternionRotationAngleDistanceDegree(q, s));
      SpzLog("[SPZ] Axis angle distance (degrees): %g", quaternionAxisAngleDistanceDegree(q, s));
#endif
    }
#endif

    if (o.fast_rot_quantization) {
      for (int j = 0; j < 3; ++j) packed.rotations[i * 3 + j] = r0[j];
      continue;
    } else {
      // Optimize the quaternion representation by exploring a neighborhood around the initial quantization
      // to find the best 8-bit representation that minimizes reconstruction error
            
      // Track the best representation found so far
      float err_min = 1;  // Start with maximum possible error (1.0 means opposite quaternions)
      int radius_err_min = -1;  // Best error found so far
      int radius_tested_max = -1;  // Max radius at which point was tested
      uint8_t r_best[3] = {r0[0], r0[1], r0[2]};  // Start with the initial quantization
      
      // Progressively search larger neighborhoods until we find a good enough representation
      // or reach the maximum radius
      for (int radius = 0; radius <= radius_max && (radius == 0 || err_min >= errorThreshold[radius - 1]); ++radius) {
#if VERBOSE >= 3
        SpzLog("[SPZ] radius=%d", radius);
#endif
        if (radius >  radius_tested_max) {
          radius_tested_max = radius;
        }
        // Special case for radius=0: just evaluate the initial quantization
        if (radius == 0) {
          Quat4f s = unpackQuaternion(r0);
          float err = quaternionReconstructionError(q, s);
          err_min = err;
          radius_err_min = radius;
#if VERBOSE >= 3
          SpzLog("[SPZ] err_min=%g", err_min);
#endif
          continue;
        }
        
        // For radius > 0, explore the surface of a cube with the given radius
        // We only need to check the surface since we've already checked all points
        // with smaller radii in previous iterations
        for (int dim = 0; dim < 3; ++dim) {
          // For each dimension, check both the positive and negative faces
          for (int val : {-radius, radius}) {
            // Set up the other two dimensions for iteration
            int nextDim = (dim + 1) % 3;
            int lastDim = (dim + 2) % 3;
            
            // Define range limits for the other two dimensions to avoid duplicate processing
            // We process faces in a specific order to avoid checking the same point multiple times:
            // - x-faces: full y and z ranges
            // - y-faces: full z range but limited x range (excluding edges already covered by x-faces)
            // - z-faces: limited x and y ranges (excluding edges already covered by x and y faces)
            int nextMin = (dim < 2) ? -radius : -radius+1;
            int nextMax = (dim < 2) ? radius : radius-1;
            int lastMin = (dim < 1) ? -radius : -radius+1;
            int lastMax = (dim < 1) ? radius : radius-1;
            
            // Iterate over all points on this face
            for (int face_i = nextMin; face_i <= nextMax; ++face_i) {
              for (int face_j = lastMin; face_j <= lastMax; ++face_j) {
                // Calculate the coordinates of this point relative to r0
                int coords[3] = {0, 0, 0};
                coords[dim] = r0[dim] + val;
                coords[nextDim] = r0[nextDim] + face_i;
                coords[lastDim] = r0[lastDim] + face_j;
                
                // Skip points outside the valid 8-bit range
                if (coords[0] < 0 || coords[0] > 255 || 
                    coords[1] < 0 || coords[1] > 255 ||
                    coords[2] < 0 || coords[2] > 255) {
                  continue;
                }
                
                // Create the candidate quaternion representation
                uint8_t r[3];
                for (int j = 0; j < 3; ++j) r[j] = static_cast<uint8_t>(coords[j]);
                
                // Unpack it to evaluate its quality
                Quat4f s = unpackQuaternion(r);
                
#if VERBOSE >= 3
                {
                  SpzLog("[SPZ] Packed");
                  SpzLog("[SPZ] %d, %d, %d", (int)r[0], (int)r[1], (int)r[2]);
                  SpzLog("[SPZ] Unpacked normalized");
                  SpzLog("[SPZ] %g, %g, %g, %g (scalar)", s[0], s[1], s[2], s[3]);

                  SpzLog("[SPZ] Distances:");
                  SpzLog("[SPZ] Reconstruction error: %g", quaternionReconstructionError(q, s));
                  SpzLog("[SPZ] Geodesic quaternion distance (degrees): %g", quaternionGeodesicDistanceDegree(q, s));
                  SpzLog("[SPZ] Rotation angle distance (degrees): %g", quaternionRotationAngleDistanceDegree(q, s));
                  SpzLog("[SPZ] Axis angle distance (degrees): %g", quaternionAxisAngleDistanceDegree(q, s));
                }
#endif
                // Calculate reconstruction error and update best solution if improved
                float err = quaternionReconstructionError(q, s);
                if (err < err_min) {
                  err_min = err;
                  radius_err_min = radius;
                  for (int j = 0; j < 3; ++j) r_best[j] = r[j];
#if VERBOSE >= 3
                  SpzLog("[SPZ] err_min=%g", err_min);
#endif
                }
              }
            }
          }
        }
      }
#if VERBOSE >= 1
      if (!o.fast_rot_quantization) {
        // Unpack it to evaluate its quality
        Quat4f s = unpackQuaternion(r_best);
        float dist = quaternionGeodesicDistanceDegree(q, s);
        averageGeodesicDistanceDegreeBest += dist;
#if VERBOSE >= 2
        SpzLog("[SPZ] Winner found at radius=%d (packed):", radius_err_min);
        SpzLog("[SPZ] %d, %d, %d", (int)r_best[0], (int)r_best[1], (int)r_best[2]);
        SpzLog("[SPZ] Distances:");
        SpzLog("[SPZ] Reconstruction error: %g", quaternionGeodesicDistanceDegreeToReconstructionError(dist));
        SpzLog("[SPZ] Geodesic quaternion distance (degrees): %g", dist);
        SpzLog("[SPZ] Rotation angle distance (degrees): %g", quaternionRotationAngleDistanceDegree(q, s));
        SpzLog("[SPZ] Axis angle distance (degrees): %g", quaternionAxisAngleDistanceDegree(q, s));
#endif
      }
#endif
      for (int j = 0; j < 3; ++j) packed.rotations[i * 3 + j] = r_best[j];
      numPointsMaxRadius[radius_tested_max] += 1;
      numPointsBestRadius[radius_err_min] += 1;
    }
  }
#if VERBOSE >= 1
  averageGeodesicDistanceDegreeFast /= numPoints;
  averageGeodesicDistanceDegreeBest /= numPoints;
  SpzLog("[SPZ] Rotation quantization statistics:");
  SpzLog("[SPZ] Average geodesic distance in degrees with fast quantization: %g", averageGeodesicDistanceDegreeFast);
  if (!o.fast_rot_quantization) {
    SpzLog("[SPZ] Average geodesic distance in degrees with optimized quantization: %g", averageGeodesicDistanceDegreeBest);
    for (int radius = 0; radius <= radius_max; ++radius) {
      SpzLog("[SPZ] Points tested up to radius %d: %.2f%% (%d)", radius, 100.0f * numPointsMaxRadius[radius] / numPoints, numPointsMaxRadius[radius]);
    }
    SpzLog("[SPZ] Best radius found:");
    for (int radius = 0; radius <= radius_max; ++radius) {
      SpzLog("[SPZ] Best radius %d: %.2f%% (%d)", radius, 100.0f * numPointsBestRadius[radius] / numPoints, numPointsBestRadius[radius]);
    }
  }
#endif

  for (size_t i = 0; i < numPoints; i++) {
    // Apply sigmoid activation to alpha
    packed.alphas[i] = toUint8(sigmoid(g.alphas[i]) * 255.0f);
  }

  for (size_t i = 0; i < numPoints * 3; i++) {
    // Convert SH DC component to wide RGB (allowing values that are a bit above 1 and below 0).
    packed.colors[i] = toUint8(g.colors[i] * (colorScale * 255.0f) + (0.5f * 255.0f));
  }

  if (g.shDegree > 0) {
    // Spherical harmonics quantization parameters. The data format uses 8 bits per coefficient, but
    // when packing, we can quantize to fewer bits for better compression.
    constexpr int32_t sh1Bits = 5;
    constexpr int32_t shRestBits = 4;
    const int32_t shPerPoint = dimForDegree(g.shDegree) * 3;
    for (size_t i = 0; i < numPoints * shPerPoint; i += shPerPoint) {
      size_t j = 0, k = 0;
      for (; j < 9; j += 3, k++) {  // There are 9 (3 * 3) coefficients for degree 1
        packed.sh[i + j + 0] = quantizeSH(c.flipSh[k] * g.sh[i + j + 0], 1 << (8 - sh1Bits));
        packed.sh[i + j + 1] = quantizeSH(c.flipSh[k] * g.sh[i + j + 1], 1 << (8 - sh1Bits));
        packed.sh[i + j + 2] = quantizeSH(c.flipSh[k] * g.sh[i + j + 2], 1 << (8 - sh1Bits));
      }
      for (; j < shPerPoint; j += 3, k++) {
        packed.sh[i + j + 0] = quantizeSH(c.flipSh[k] * g.sh[i + j + 0], 1 << (8 - shRestBits));
        packed.sh[i + j + 1] = quantizeSH(c.flipSh[k] * g.sh[i + j + 1], 1 << (8 - shRestBits));
        packed.sh[i + j + 2] = quantizeSH(c.flipSh[k] * g.sh[i + j + 2], 1 << (8 - shRestBits));
      }
    }
  }

  return packed;
}

UnpackedGaussian PackedGaussian::unpack(
  bool usesFloat16, int32_t fractionalBits, const CoordinateConverter &c) const {
  UnpackedGaussian result;
  if (usesFloat16) {
    // Decode legacy float16 format. We can remove this at some point as it was never released.
    const auto *halfData = reinterpret_cast<const Half *>(position.data());
    for (size_t i = 0; i < 3; i++) {
      result.position[i] = c.flipP[i] * halfToFloat(halfData[i]);
    }
  } else {
    // Decode 24-bit fixed point coordinates
    float scale = 1.0 / (1 << fractionalBits);
    for (size_t i = 0; i < 3; i++) {
      int32_t fixed32 = position[i * 3 + 0];
      fixed32 |= position[i * 3 + 1] << 8;
      fixed32 |= position[i * 3 + 2] << 16;
      fixed32 |= (fixed32 & 0x800000) ? 0xff000000 : 0;  // sign extension
      result.position[i] = c.flipP[i] * static_cast<float>(fixed32) * scale;
    }
  }

  for (size_t i = 0; i < 3; i++) {
    result.scale[i] = (scale[i] / 16.0f - 10.0f);
  }

  const uint8_t *r = &rotation[0];
  Vec3f xyz = times(
    plus(
      times(
        Vec3f{static_cast<float>(r[0]), static_cast<float>(r[1]), static_cast<float>(r[2])},
        1.0f / 127.5f),
      Vec3f{-1, -1, -1}),
    c.flipQ);
  std::copy(xyz.data(), xyz.data() + 3, &result.rotation[0]);
  // Compute the real component - we know the quaternion is normalized and w is non-negative
  result.rotation[3] = std::sqrt(std::max(0.0f, 1.0f - squaredNorm(xyz)));

  result.alpha = invSigmoid(alpha / 255.0f);

  for (size_t i = 0; i < 3; i++) {
    result.color[i] = ((color[i] / 255.0f) - 0.5f) / colorScale;
  }

  for (size_t i = 0; i < 15; i++) {
    result.shR[i] = c.flipSh[i] * unquantizeSH(shR[i]);
    result.shG[i] = c.flipSh[i] * unquantizeSH(shG[i]);
    result.shB[i] = c.flipSh[i] * unquantizeSH(shB[i]);
  }

  return result;
}

PackedGaussian PackedGaussians::at(int32_t i) const {
  PackedGaussian result;
  int32_t positionBits = usesFloat16() ? 6 : 9;
  int32_t start3 = i * 3;
  const auto *p = &positions[i * positionBits];
  std::copy(p, p + positionBits, result.position.data());
  std::copy(&scales[start3], &scales[start3] + 3, result.scale.data());
  std::copy(&rotations[start3], &rotations[start3] + 3, result.rotation.data());
  std::copy(&colors[start3], &colors[start3] + 3, result.color.data());
  result.alpha = alphas[i];

  int32_t shDim = dimForDegree(shDegree);
  const auto *sh = &this->sh[i * shDim * 3];
  for (int32_t j = 0; j < shDim; ++j, sh += 3) {
    result.shR[j] = sh[0];
    result.shG[j] = sh[1];
    result.shB[j] = sh[2];
  }
  for (int32_t j = shDim; j < 15; ++j) {
    result.shR[j] = 128;
    result.shG[j] = 128;
    result.shB[j] = 128;
  }

  return result;
}

UnpackedGaussian PackedGaussians::unpack(int32_t i, const CoordinateConverter &c) const {
  return at(i).unpack(usesFloat16(), fractionalBits, c);
}

bool PackedGaussians::usesFloat16() const { return positions.size() == numPoints * 3 * 2; }

GaussianCloud unpackGaussians(const PackedGaussians &packed, const UnpackOptions &o) {
  const int32_t numPoints = packed.numPoints;
  const int32_t shDim = dimForDegree(packed.shDegree);
  const bool usesFloat16 = packed.usesFloat16();
  if (!checkSizes(packed, numPoints, shDim, usesFloat16)) {
    return {};
  }

  GaussianCloud result;
  result.numPoints = packed.numPoints;
  result.shDegree = packed.shDegree;
  result.antialiased = packed.antialiased;
  result.positions.resize(numPoints * 3);
  result.scales.resize(numPoints * 3);
  result.rotations.resize(numPoints * 4);
  result.alphas.resize(numPoints);
  result.colors.resize(numPoints * 3);
  result.sh.resize(numPoints * shDim * 3);

  if (usesFloat16) {
    // Decode legacy float16 format. We can remove this at some point as it was never released.
    const auto *halfData = reinterpret_cast<const Half *>(packed.positions.data());
    for (size_t i = 0; i < numPoints * 3; i++) {
      result.positions[i] = halfToFloat(halfData[i]);
    }
  } else {
    // Decode 24-bit fixed point coordinates
    float scale = 1.0 / (1 << packed.fractionalBits);
    for (size_t i = 0; i < numPoints * 3; i++) {
      int32_t fixed32 = packed.positions[i * 3 + 0];
      fixed32 |= packed.positions[i * 3 + 1] << 8;
      fixed32 |= packed.positions[i * 3 + 2] << 16;
      fixed32 |= (fixed32 & 0x800000) ? 0xff000000 : 0;  // sign extension
      result.positions[i] = static_cast<float>(fixed32) * scale;
    }
  }

  for (size_t i = 0; i < numPoints * 3; i++) {
    result.scales[i] = packed.scales[i] / 16.0f - 10.0f;
  }

  for (size_t i = 0; i < numPoints; i++) {
    const uint8_t *r = &packed.rotations[i * 3];
    Vec3f xyz = plus(
      times(
        Vec3f{static_cast<float>(r[0]), static_cast<float>(r[1]), static_cast<float>(r[2])},
        1.0f / 127.5f),
      Vec3f{-1, -1, -1});
    std::copy(xyz.data(), xyz.data() + 3, &result.rotations[i * 4]);
    // Compute the real component - we know the quaternion is normalized and w is non-negative
    result.rotations[i * 4 + 3] = std::sqrt(std::max(0.0f, 1.0f - squaredNorm(xyz)));
  }

  for (size_t i = 0; i < numPoints; i++) {
    result.alphas[i] = invSigmoid(packed.alphas[i] / 255.0f);
  }

  for (size_t i = 0; i < numPoints * 3; i++) {
    result.colors[i] = ((packed.colors[i] / 255.0f) - 0.5f) / colorScale;
  }

  for (size_t i = 0; i < packed.sh.size(); i++) {
    result.sh[i] = unquantizeSH(packed.sh[i]);
  }

  result.convertCoordinates(CoordinateSystem::RUB, o.to);
  return result;
}

void serializePackedGaussians(const PackedGaussians &packed, std::ostream *out) {
  PackedGaussiansHeader header;
  header.numPoints = static_cast<uint32_t>(packed.numPoints);
  header.shDegree = static_cast<uint8_t>(packed.shDegree);
  header.fractionalBits = static_cast<uint8_t>(packed.fractionalBits);
  header.flags = static_cast<uint8_t>(packed.antialiased ? FlagAntialiased : 0);
  out->write(reinterpret_cast<const char *>(&header), sizeof(header));
  out->write(reinterpret_cast<const char *>(packed.positions.data()), countBytes(packed.positions));
  out->write(reinterpret_cast<const char *>(packed.alphas.data()), countBytes(packed.alphas));
  out->write(reinterpret_cast<const char *>(packed.colors.data()), countBytes(packed.colors));
  out->write(reinterpret_cast<const char *>(packed.scales.data()), countBytes(packed.scales));
  out->write(reinterpret_cast<const char *>(packed.rotations.data()), countBytes(packed.rotations));
  out->write(reinterpret_cast<const char *>(packed.sh.data()), countBytes(packed.sh));
}

PackedGaussians deserializePackedGaussians(std::istream &in) {
  constexpr int32_t maxPointsToRead = 10000000;

  PackedGaussiansHeader header;
  in.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (!in || header.magic != PackedGaussiansHeader().magic) {
    SpzLog("[SPZ ERROR] deserializePackedGaussians: header not found");
    return {};
  }
  if (header.version < 1 || header.version > 2) {
    SpzLog("[SPZ ERROR] deserializePackedGaussians: version not supported: %d", header.version);
    return {};
  }
  if (header.numPoints > maxPointsToRead) {
    SpzLog("[SPZ ERROR] deserializePackedGaussians: Too many points: %d", header.numPoints);
    return {};
  }
  if (header.shDegree > 3) {
    SpzLog("[SPZ ERROR] deserializePackedGaussians: Unsupported SH degree: %d", header.shDegree);
    return {};
  }
  const int32_t numPoints = header.numPoints;
  const int32_t shDim = dimForDegree(header.shDegree);
  const bool usesFloat16 = header.version == 1;
  PackedGaussians result;
  result.numPoints = numPoints;
  result.shDegree = header.shDegree;
  result.fractionalBits = header.fractionalBits;
  result.antialiased = (header.flags & FlagAntialiased) != 0;
  result.positions.resize(numPoints * 3 * (usesFloat16 ? 2 : 3));
  result.scales.resize(numPoints * 3);
  result.rotations.resize(numPoints * 3);
  result.alphas.resize(numPoints);
  result.colors.resize(numPoints * 3);
  result.sh.resize(numPoints * shDim * 3);
  in.read(reinterpret_cast<char *>(result.positions.data()), countBytes(result.positions));
  in.read(reinterpret_cast<char *>(result.alphas.data()), countBytes(result.alphas));
  in.read(reinterpret_cast<char *>(result.colors.data()), countBytes(result.colors));
  in.read(reinterpret_cast<char *>(result.scales.data()), countBytes(result.scales));
  in.read(reinterpret_cast<char *>(result.rotations.data()), countBytes(result.rotations));
  in.read(reinterpret_cast<char *>(result.sh.data()), countBytes(result.sh));
  if (!in) {
    SpzLog("[SPZ ERROR] deserializePackedGaussians: read error");
    return {};
  }
  return result;
}

bool saveSpz(const GaussianCloud &g, const PackOptions &o, std::vector<uint8_t> *out) {
  std::string data;
  {
    PackedGaussians packed = packGaussians(g, o);
    std::stringstream ss;
    serializePackedGaussians(packed, &ss);
    data = ss.str();
  }
  return compressGzipped(reinterpret_cast<const uint8_t *>(data.data()), data.size(), out);
}

PackedGaussians loadSpzPacked(const uint8_t *data, int32_t size) {
  std::string decompressed;
  if (!decompressGzipped(data, size, &decompressed))
    return {};
  std::stringstream stream(std::move(decompressed));
  return deserializePackedGaussians(stream);
}

PackedGaussians loadSpzPacked(const std::vector<uint8_t> &data) {
  return loadSpzPacked(data.data(), static_cast<int>(data.size()));
}

PackedGaussians loadSpzPacked(const std::string &filename) {
  std::ifstream in(filename, std::ios::binary | std::ios::ate);
  if (!in.good())
    return {};
  std::vector<uint8_t> data(in.tellg());
  in.seekg(0, std::ios::beg);
  in.read(reinterpret_cast<char *>(data.data()), data.size());
  if (!in.good()) {
    return {};
  }
  return loadSpzPacked(data);
}

GaussianCloud loadSpz(const std::vector<uint8_t> &data, const UnpackOptions &o) {
  return unpackGaussians(loadSpzPacked(data), o);
}

bool saveSpz(const GaussianCloud &g, const PackOptions &o, const std::string &filename) {
  std::vector<uint8_t> data;
  if (!saveSpz(g, o, &data)) {
    return false;
  }
  std::ofstream out(filename, std::ios::binary | std::ios::out);
  out.write(reinterpret_cast<const char *>(data.data()), data.size());
  out.close();
  return out.good();
}

GaussianCloud loadSpz(const std::string &filename, const UnpackOptions &o) {
  std::ifstream in(filename, std::ios::binary | std::ios::ate);
  if (!in.good()) {
    SpzLog("[SPZ ERROR] Unable to open: %s", filename.c_str());
    return {};
  }
  std::vector<uint8_t> data(in.tellg());
  in.seekg(0, std::ios::beg);
  in.read(reinterpret_cast<char *>(data.data()), data.size());
  in.close();
  if (!in.good()) {
    SpzLog("[SPZ ERROR] Unable to load data from: %s", filename.c_str());
    return {};
  }
  return loadSpz(data, o);
}

// Helper function to get a line from the file, skipping comment lines
bool getLineSkipComments(std::ifstream& in, std::string& line) {
  while (std::getline(in, line)) {
    if (line.compare(0, 8, "comment ") != 0) {
      return true;
    }
  }
  return false;
}

GaussianCloud loadSplatFromPly(const std::string &filename, const UnpackOptions &o) {
  SpzLog("[SPZ] Loading: %s", filename.c_str());
  std::ifstream in(filename, std::ios::binary);
  if (!in.good()) {
    SpzLog("[SPZ ERROR] Unable to open: %s", filename.c_str());
    in.close();
    return {};
  }
  std::string line;
  std::getline(in, line);
  if (line != "ply") {
    SpzLog("[SPZ ERROR] %s: not a .ply file", filename.c_str());
    in.close();
    return {};
  }
  std::getline(in, line);
  if (line != "format binary_little_endian 1.0") {
    SpzLog("[SPZ ERROR] %s: unsupported .ply format", filename.c_str());
    in.close();
    return {};
  }
  if (!getLineSkipComments(in, line) || line.find("element vertex ") != 0) {
    SpzLog("[SPZ ERROR] %s: missing vertex count", filename.c_str());
    in.close();
    return {};
  }
  int32_t numPoints = std::stoi(line.substr(std::strlen("element vertex ")));
  if (numPoints <= 0 || numPoints > 10 * 1024 * 1024) {
    SpzLog("[SPZ ERROR] %s: invalid vertex count: %d", filename.c_str(), numPoints);
    in.close();
    return {};
  }

  SpzLog("[SPZ] Loading %d points", numPoints);
  std::unordered_map<std::string, int> fields;  // name -> index
  for (int32_t i = 0;; i++) {
    if (!getLineSkipComments(in, line)) {
      SpzLog("[SPZ ERROR] %s: unexpected end of file", filename.c_str());
      in.close();
      return {};
    }
    
    if (line == "end_header")
      break;

    if (line.find("property float ") != 0) {
      SpzLog("[SPZ ERROR] %s: unsupported property data type: %s", filename.c_str(), line.c_str());
      in.close();
      return {};
    }
    std::string name = line.substr(std::strlen("property float "));
    fields[name] = i;
  }

  // Returns the index for a given field name, ensuring the name exists.
  const auto index = [&fields](const std::string &name) {
    const auto &itr = fields.find(name);
    if (itr == fields.end()) {
      SpzLog("[SPZ ERROR] Missing field: %s", name.c_str());
      return -1;
    }
    return itr->second;
  };

  const std::vector<int> positionIdx = {index("x"), index("y"), index("z")};
  const std::vector<int> scaleIdx = {index("scale_0"), index("scale_1"), index("scale_2")};
  const std::vector<int> rotIdx = {index("rot_1"), index("rot_2"), index("rot_3"), index("rot_0")};
  const std::vector<int> alphaIdx = {index("opacity")};
  const std::vector<int> colorIdx = {index("f_dc_0"), index("f_dc_1"), index("f_dc_2")};

  // Check that only valid indices were returned.
  for (auto idx : positionIdx) {
    if (idx < 0) {
      in.close();
      return {};
    }
  }
  for (auto idx : scaleIdx) {
    if (idx < 0) {
      in.close();
      return {};
    }
  }
  for (auto idx : rotIdx) {
    if (idx < 0) {
      in.close();
      return {};
    }
  }
  for (auto idx : alphaIdx) {
    if (idx < 0) {
      in.close();
      return {};
    }
  }
  for (auto idx : colorIdx) {
    if (idx < 0) {
      in.close();
      return {};
    }
  }

  // Spherical harmonics are optional and variable in size (depending on degree)
  std::vector<int> shIdx;
  for (int32_t i = 0; i < 45; i++) {
    const auto &itr = fields.find("f_rest_" + std::to_string(i));
    if (itr == fields.end())
      break;
    shIdx.push_back(itr->second);
  }
  const int32_t shDim = static_cast<int>(shIdx.size() / 3);

  std::vector<float> values(numPoints * fields.size());
  in.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(float));
  if (!in.good()) {
    SpzLog("[SPZ ERROR] Unable to load data from: %s", filename.c_str());
    in.close();
    return {};
  }
  in.close();

  GaussianCloud result;
  result.numPoints = numPoints;
  result.shDegree = degreeForDim(shDim);
  result.positions.reserve(numPoints * 3);
  result.scales.reserve(numPoints * 3);
  result.rotations.reserve(numPoints * 4);
  result.alphas.reserve(numPoints * 1);
  result.colors.reserve(numPoints * 3);
  for (size_t i = 0; i < values.size(); i += fields.size()) {
    for (int32_t j = 0; j < positionIdx.size(); j++) {
      result.positions.push_back(values[i + positionIdx[j]]);
    }
    for (int32_t j = 0; j < scaleIdx.size(); j++) {
      result.scales.push_back(values[i + scaleIdx[j]]);
    }
    for (int32_t j = 0; j < rotIdx.size(); j++) {
      result.rotations.push_back(values[i + rotIdx[j]]);
    }
    for (int32_t j = 0; j < alphaIdx.size(); j++) {
      result.alphas.push_back(values[i + alphaIdx[j]]);
    }
    for (int32_t j = 0; j < colorIdx.size(); j++) {
      result.colors.push_back(values[i + colorIdx[j]]);
    }
    // Convert from [N,C,S] to [N,S,C] (where C is color channel, S is SH coeff).
    for (int32_t j = 0; j < shDim; j++) {
      result.sh.push_back(values[i + shIdx[j]]);
      result.sh.push_back(values[i + shIdx[j + shDim]]);
      result.sh.push_back(values[i + shIdx[j + 2 * shDim]]);
    }
  }

  result.convertCoordinates(CoordinateSystem::RDF, o.to);
  return result;
}

bool saveSplatToPly(const GaussianCloud &data, const PackOptions &o, const std::string &filename) {
  const int32_t N = data.numPoints;
  CHECK_EQ(data.positions.size(), N * 3);
  CHECK_EQ(data.scales.size(), N * 3);
  CHECK_EQ(data.rotations.size(), N * 4);
  CHECK_EQ(data.alphas.size(), N);
  CHECK_EQ(data.colors.size(), N * 3);
  const int32_t shDim = static_cast<int>(data.sh.size() / N / 3);
  const int32_t D = 17 + shDim * 3;

  CoordinateConverter c = coordinateConverter(o.from, CoordinateSystem::RDF);

  std::vector<float> values(N * D, 0.0f);
  int32_t outIdx = 0, i3 = 0, i4 = 0;
  for (int32_t i = 0; i < N; i++) {
    // Position (x, y, z)
    values[outIdx++] = c.flipP[0] * data.positions[i3 + 0];
    values[outIdx++] = c.flipP[1] * data.positions[i3 + 1];
    values[outIdx++] = c.flipP[2] * data.positions[i3 + 2];
    // Normals (nx, ny, nz): these are always zero, but some viewers expect them to be present
    outIdx += 3;
    // Color (r, g, b): DC component for spherical harmonics
    values[outIdx++] = data.colors[i3 + 0];
    values[outIdx++] = data.colors[i3 + 1];
    values[outIdx++] = data.colors[i3 + 2];
    // Spherical harmonics: Interleave so the coefficients are the fastest-changing axis and
    // the channel (r, g, b) is slower-changing axis.
    for (int32_t j = 0; j < shDim; j++) {
      values[outIdx++] = c.flipSh[j] * data.sh[(i * shDim + j) * 3];
    }
    for (int32_t j = 0; j < shDim; j++) {
      values[outIdx++] = c.flipSh[j] * data.sh[(i * shDim + j) * 3 + 1];
    }
    for (int32_t j = 0; j < shDim; j++) {
      values[outIdx++] = c.flipSh[j] * data.sh[(i * shDim + j) * 3 + 2];
    }
    // Alpha
    values[outIdx++] = data.alphas[i];
    // Scale (sx, sy, sz)
    values[outIdx++] = data.scales[i3 + 0];
    values[outIdx++] = data.scales[i3 + 1];
    values[outIdx++] = data.scales[i3 + 2];
    // Rotation (qw, qx, qy, qz)
    values[outIdx++] = data.rotations[i4 + 3];
    values[outIdx++] = c.flipQ[0] * data.rotations[i4 + 0];
    values[outIdx++] = c.flipQ[1] * data.rotations[i4 + 1];
    values[outIdx++] = c.flipQ[2] * data.rotations[i4 + 2];
    i3 += 3;
    i4 += 4;
  }
  CHECK_EQ(outIdx, values.size());

  std::ofstream out(filename, std::ios::binary);
  if (!out.good()) {
    SpzLog("[SPZ ERROR] Unable to open for writing: %s", filename.c_str());
    return false;
  }
  out << "ply\n";
  out << "format binary_little_endian 1.0\n";
  out << "element vertex " << N << "\n";
  out << "property float x\n";
  out << "property float y\n";
  out << "property float z\n";
  out << "property float nx\n";
  out << "property float ny\n";
  out << "property float nz\n";
  out << "property float f_dc_0\n";
  out << "property float f_dc_1\n";
  out << "property float f_dc_2\n";
  for (int32_t i = 0; i < shDim * 3; i++) {
    out << "property float f_rest_" << i << "\n";
  }
  out << "property float opacity\n";
  out << "property float scale_0\n";
  out << "property float scale_1\n";
  out << "property float scale_2\n";
  out << "property float rot_0\n";
  out << "property float rot_1\n";
  out << "property float rot_2\n";
  out << "property float rot_3\n";
  out << "end_header\n";
  out.write(reinterpret_cast<char *>(values.data()), values.size() * sizeof(float));
  out.close();
  if (!out.good()) {
    SpzLog("[SPZ ERROR] Failed to write to: %s", filename.c_str());
    return false;
  }
  return true;
}

}  // namespace spz
