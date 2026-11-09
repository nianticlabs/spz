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


#include "temporal-niantic.h"
#include "extension-frames.h"
#include "load-spz.h"
#include "quantized-planes.h"
#include "splat-utils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>

namespace spz {
namespace {
constexpr size_t kPositionComponents = 3;
constexpr size_t kQuaternionComponents = 4;
// Window channels hold [center, ln sigma, peak opacity], optionally followed by an interval [lo, hi).
constexpr size_t kWindowComponents = 3;
constexpr size_t kWindowComponentsWithInterval = 5;
constexpr size_t kWindowCenter = 0;
constexpr size_t kWindowLogSigma = 1;
constexpr size_t kWindowPeak = 2;
constexpr size_t kWindowLo = 3;
constexpr size_t kWindowHi = 4;
constexpr uint8_t kMaxPolynomialDegree = 3;
constexpr size_t kMaxSkinInfluences = 8;
// Skin channels hold an index and a weight per influence.
constexpr size_t kSkinValuesPerInfluence = 2;
// Driver keyframes hold translation xyz, then the rotation quaternion from this offset.
constexpr size_t kDriverRotation = 3;
// A ping-pong cycle plays the clip forward, then backward.
constexpr double kPingPongPasses = 2.0;
constexpr uint8_t kMaxLoop = static_cast<uint8_t>(SpzTemporalLoop::PingPong);
constexpr uint8_t kMaxEncoding = static_cast<uint8_t>(SpzChannelEncoding::Step);
constexpr size_t kMaxNameBytes = std::numeric_limits<uint16_t>::max();
constexpr char kFirstPrintable = 0x20;
constexpr char kLastPrintable = 0x7E;
// Records are keyed so a file holds at most one; the key itself carries no information.
constexpr char kRecordKey[] = "temporal";
// Time bins writers use when a clip has no nominal frame rate, and the most they ever use.
constexpr size_t kDefaultTimeBins = 64;
constexpr size_t kMaxTimeBins = 4096;
// Alphas are stored as logits, which can't represent 0 or 1 exactly.
constexpr float kMinAlpha = 1e-6f;
// A window this many clip durations wide changes opacity by under 1e-6 within the clip, so the
// builder treats it as always visible.
constexpr double kAlwaysVisibleSigmaRatio = 1000.0;
// Above this cosine, slerp falls back to normalized lerp, which is then accurate and stable.
constexpr float kSlerpLerpThreshold = 0.9995f;

using Vec3 = std::array<float, 3>;
using Quat = std::array<float, 4>;  // x, y, z, w — the order clouds store rotations in

// ---- Quaternion helpers (xyzw) ----

float dot4(const Quat& a, const Quat& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

// Returns q / |q|, or q unchanged if it is zero.
Quat normalize4(const Quat& q) {
  const float n = std::sqrt(dot4(q, q));
  if (n == 0.0f) return q;
  return {q[0] / n, q[1] / n, q[2] / n, q[3] / n};
}

Quat quatAt(const float* xyzw) { return {xyzw[0], xyzw[1], xyzw[2], xyzw[3]}; }

// Hamilton product a ⊗ b.
Quat multiply(const Quat& a, const Quat& b) {
  return {a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
          a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
          a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
          a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]};
}

// Rotates v by the unit quaternion q.
Vec3 rotate(const Quat& q, const Vec3& v) {
  // v + 2w (u × v) + 2 u × (u × v), with u the vector part of q.
  const Vec3 u = {q[0], q[1], q[2]};
  const Vec3 uv = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  const Vec3 uuv = {u[1] * uv[2] - u[2] * uv[1], u[2] * uv[0] - u[0] * uv[2],
                    u[0] * uv[1] - u[1] * uv[0]};
  return {v[0] + 2.0f * (q[3] * uv[0] + uuv[0]), v[1] + 2.0f * (q[3] * uv[1] + uuv[1]),
          v[2] + 2.0f * (q[3] * uv[2] + uuv[2])};
}

// Shortest-path spherical interpolation between unit quaternions.
Quat slerp(const Quat& a, Quat b, float u) {
  float cosine = dot4(a, b);
  if (cosine < 0.0f) {
    b = {-b[0], -b[1], -b[2], -b[3]};
    cosine = -cosine;
  }
  if (cosine > kSlerpLerpThreshold) {
    return normalize4({a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1]),
                       a[2] + u * (b[2] - a[2]), a[3] + u * (b[3] - a[3])});
  }
  const float angle = std::acos(cosine);
  const float sinAngle = std::sin(angle);
  const float wa = std::sin((1.0f - u) * angle) / sinAngle;
  const float wb = std::sin(u * angle) / sinAngle;
  return {wa * a[0] + wb * b[0], wa * a[1] + wb * b[1], wa * a[2] + wb * b[2],
          wa * a[3] + wb * b[3]};
}

using Mat3 = std::array<std::array<float, 3>, 3>;  // [row][column]

Mat3 toMatrix(const Quat& q) {
  const float x = q[0], y = q[1], z = q[2], w = q[3];
  return {{{1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)},
           {2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)},
           {2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)}}};
}

// Shepperd's method: picks the numerically largest of w, x, y, z to divide by.
Quat fromMatrix(const Mat3& m) {
  const float trace = m[0][0] + m[1][1] + m[2][2];
  Quat q;
  if (trace > 0.0f) {
    const float s = std::sqrt(trace + 1.0f) * 2.0f;
    q = {(m[2][1] - m[1][2]) / s, (m[0][2] - m[2][0]) / s, (m[1][0] - m[0][1]) / s, 0.25f * s};
  } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
    const float s = std::sqrt(1.0f + m[0][0] - m[1][1] - m[2][2]) * 2.0f;
    q = {0.25f * s, (m[0][1] + m[1][0]) / s, (m[0][2] + m[2][0]) / s, (m[2][1] - m[1][2]) / s};
  } else if (m[1][1] > m[2][2]) {
    const float s = std::sqrt(1.0f + m[1][1] - m[0][0] - m[2][2]) * 2.0f;
    q = {(m[0][1] + m[1][0]) / s, 0.25f * s, (m[1][2] + m[2][1]) / s, (m[0][2] - m[2][0]) / s};
  } else {
    const float s = std::sqrt(1.0f + m[2][2] - m[0][0] - m[1][1]) * 2.0f;
    q = {(m[0][2] + m[2][0]) / s, (m[1][2] + m[2][1]) / s, 0.25f * s, (m[1][0] - m[0][1]) / s};
  }
  return normalize4(q);
}

Mat3 multiply(const Mat3& a, const Mat3& b) {
  Mat3 r{};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k) r[i][j] += a[i][k] * b[k][j];
  return r;
}

Mat3 transposed(const Mat3& m) {
  Mat3 r{};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j) r[i][j] = m[j][i];
  return r;
}

float sigmoid(float x) { return 1.0f / (1.0f + std::exp(-x)); }

float alphaLogit(double alpha) {
  const double clamped = std::clamp(alpha, double{kMinAlpha}, 1.0 - kMinAlpha);
  return static_cast<float>(std::log(clamped / (1.0 - clamped)));
}

void logStatic(const std::string& problem) {
  SpzLog("[SPZ ERROR] Temporal record: %s; treating the cloud as static", problem.c_str());
}

// ---- Channels ----

bool isRequired(const SpzTemporalChannel& channel) {
  return (static_cast<uint8_t>(channel.flags) &
          static_cast<uint8_t>(SpzTemporalChannelFlags::Required)) != 0;
}

bool isSupported(const SpzTemporalChannel& channel) {
  using T = SpzTemporalTarget;
  using Y = SpzTemporalChannelType;
  return (channel.target == T::Position && channel.type == Y::Polynomial) ||
         (channel.target == T::Rotation && channel.type == Y::Polynomial) ||
         (channel.target == T::Opacity && channel.type == Y::Window) ||
         (channel.target == T::Transform && channel.type == Y::Skin);
}

bool isRotationPolynomial(const SpzTemporalChannel& c) {
  return c.target == SpzTemporalTarget::Rotation && c.type == SpzTemporalChannelType::Polynomial;
}

bool isPositionPolynomial(const SpzTemporalChannel& c) {
  return c.target == SpzTemporalTarget::Position && c.type == SpzTemporalChannelType::Polynomial;
}

bool isRange(const SpzTemporalChannel& c) { return c.encoding == SpzChannelEncoding::Range; }
bool isStep(const SpzTemporalChannel& c) { return c.encoding == SpzChannelEncoding::Step; }
bool isKnownEncoding(uint8_t encoding) { return encoding <= kMaxEncoding; }

// Returns why channel's encoding or its parameters are invalid, or an empty string if they are valid.
std::string encodingProblem(const SpzTemporalChannel& c) {
  switch (c.encoding) {
    case SpzChannelEncoding::Float32:
      return {};
    case SpzChannelEncoding::Range:
      if (!isQuantizedBits(c.bits)) return "range channels need bits 1-16";
      return {};
    case SpzChannelEncoding::Step:
      if (c.steps.size() != c.dim || !validSteps(c.steps))
        return "step channels need dim finite, positive steps";
      return {};
  }
  return "unknown channel encoding";
}

// Bytes a channel's data takes for points covered points, or std::nullopt if it can't be addressed
// on this platform. Precondition: the encoding is known.
std::optional<size_t> channelBytes(const SpzTemporalChannel& c, size_t points) {
  return componentMajorBytes(points, c.dim, isRange(c) ? c.bits : kFloat32Bits);
}

bool isPrintable(const std::string& s) {
  return std::all_of(s.begin(), s.end(),
                     [](char c) { return c >= kFirstPrintable && c <= kLastPrintable; });
}

// Returns why a supported channel's dim doesn't fit its layout, or an empty string if it does.
std::string layoutProblem(const SpzTemporalChannel& c) {
  switch (c.target) {
    case SpzTemporalTarget::Position:
      if (c.dim % kPositionComponents != 0 || c.dim / kPositionComponents > kMaxPolynomialDegree)
        return "position channels need dim 3, 6 or 9";
      return {};
    case SpzTemporalTarget::Rotation:
      if (c.dim % kQuaternionComponents != 0 ||
          c.dim / kQuaternionComponents > kMaxPolynomialDegree) {
        return "rotation channels need dim 4, 8 or 12";
      }
      return {};
    case SpzTemporalTarget::Opacity:
      if (c.dim != kWindowComponents && c.dim != kWindowComponentsWithInterval)
        return "opacity channels need dim 3, or 5 with an interval";
      return {};
    case SpzTemporalTarget::Transform:
      if (c.dim % kSkinValuesPerInfluence != 0 ||
          c.dim / kSkinValuesPerInfluence > kMaxSkinInfluences) {
        return "skin channels need dim 2, 4, ..., 16";
      }
      return {};
  }
  return "unsupported channel";
}

// The number of points channel covers. Precondition: its point range fits numPoints.
size_t coveredPoints(const SpzTemporalChannel& c, int32_t numPoints) {
  return c.pointRange ? c.pointRange->end - c.pointRange->begin : static_cast<size_t>(numPoints);
}

std::pair<size_t, size_t> coveredRange(const SpzTemporalChannel& c, int32_t numPoints) {
  if (c.pointRange) return {c.pointRange->begin, c.pointRange->end};
  return {0, static_cast<size_t>(numPoints)};
}

const char* rangeProblem(const SpzTemporalChannel& c, int32_t numPoints) {
  if (!c.pointRange) return nullptr;
  if (c.pointRange->begin >= c.pointRange->end) return "a channel's point range must have begin < end";
  if (numPoints < 0 || c.pointRange->end > static_cast<uint32_t>(numPoints))
    return "a channel's point range exceeds the cloud's points";
  return nullptr;
}

// Conversions permute and mix a vector's components, so a Step channel's per-component steps no
// longer line up with them. Gives every component of each vector its vector's smallest step, which
// never loosens the error bound.
void unifyStepsPerVector(size_t components, SpzTemporalChannel* c) {
  if (!isStep(*c) || c->steps.size() != c->dim) return;
  for (size_t start = 0; start + components <= c->steps.size(); start += components) {
    const auto first = c->steps.begin() + static_cast<std::ptrdiff_t>(start);
    const float smallest = *std::min_element(first, first + static_cast<std::ptrdiff_t>(components));
    std::fill(first, first + static_cast<std::ptrdiff_t>(components), smallest);
  }
}

// ---- Validation ----

// Returns why the record's fields, apart from per-point and driver data, are invalid, or an empty
// string if they are valid. Parsing checks these before decoding the record's extension frames.
std::string fieldsProblem(const SpzExtensionTemporalNiantic& r) {
  if (!std::isfinite(r.clipStart) || !std::isfinite(r.clipEnd) || !std::isfinite(r.posterTime))
    return "clip times must be finite";
  if (!(r.clipStart < r.clipEnd)) return "clipStart must be before clipEnd";
  if (r.posterTime < r.clipStart || r.posterTime >= r.clipEnd)
    return "posterTime must be in [clipStart, clipEnd)";
  if (!(r.opacityCutoff >= 0.0f && r.opacityCutoff < 1.0f)) return "opacityCutoff must be in [0, 1)";
  if (!std::isfinite(r.nominalFps) || r.nominalFps < 0.0f) return "nominalFps must be >= 0";
  if (static_cast<uint8_t>(r.loop) > kMaxLoop) return "unknown loop mode";
  if (r.channels.size() > std::numeric_limits<uint8_t>::max()) return "too many channels";
  for (size_t i = 0; i < r.channels.size(); ++i) {
    const SpzTemporalChannel& c = r.channels[i];
    if (c.dim == 0) return "channels need dim >= 1";
    std::string encodingIssue = encodingProblem(c);
    if (!encodingIssue.empty()) return encodingIssue;
    if (isSupported(c)) {
      std::string problem = layoutProblem(c);
      if (!problem.empty()) return problem;
    } else if (isRequired(c)) {
      return "required channel (target " + std::to_string(static_cast<int>(c.target)) +
             ", type " + std::to_string(static_cast<int>(c.type)) + ") isn't supported";
    }
    for (size_t j = 0; j < i; ++j) {
      if (r.channels[j].target == c.target) return "two channels animate the same target";
    }
  }
  const size_t keyframes = r.driverKeyframeCount();
  if ((r.driverCount == 0) != (keyframes == 0))
    return "drivers and driver keyframes must both be present or both absent";
  for (size_t k = 0; k < keyframes; ++k) {
    if (!std::isfinite(r.driverKeyframeTimes[k])) return "driver keyframe times must be finite";
    if (k > 0 && !(r.driverKeyframeTimes[k - 1] < r.driverKeyframeTimes[k]))
      return "driver keyframe times must be strictly increasing";
  }
  const bool skinned = std::any_of(r.channels.begin(), r.channels.end(), [](const auto& c) {
    return c.target == SpzTemporalTarget::Transform && isSupported(c);
  });
  if (skinned && r.driverCount == 0) return "a Transform channel needs drivers";
  if (!r.driverNames.empty() && r.driverNames.size() != r.driverCount)
    return "driverNames must be empty or name every driver";
  for (const std::string& name : r.driverNames) {
    if (name.empty() || name.size() > kMaxNameBytes || !isPrintable(name))
      return "driver names must be 1-65535 bytes of printable ASCII";
  }
  return {};
}

// Returns why driverTransforms doesn't fit the record's drivers, or an empty string if it does.
std::string driverTransformsProblem(const SpzExtensionTemporalNiantic& r) {
  const uint64_t expected =
    uint64_t{r.driverCount} * r.driverKeyframeCount() * kDriverTransformFloats;
  if (r.driverTransforms.size() != expected)
    return "driverTransforms must hold driverCount * keyframes * 7 floats";
  for (size_t i = 0; i < r.driverTransforms.size(); i += kDriverTransformFloats) {
    const float* f = &r.driverTransforms[i];
    if (!std::all_of(f, f + kDriverTransformFloats, [](float v) { return std::isfinite(v); }))
      return "driver transforms must be finite";
    const Quat q = quatAt(f + kDriverRotation);
    if (dot4(q, q) == 0.0f) return "driver rotations must be non-zero quaternions";
  }
  return {};
}

// Returns why the channels' per-point data doesn't fit a cloud of numPoints points, or an empty
// string if it does.
std::string channelDataProblem(const SpzExtensionTemporalNiantic& r, int32_t numPoints) {
  for (const SpzTemporalChannel& c : r.channels) {
    if (const char* problem = rangeProblem(c, numPoints)) return problem;
    if (c.values.size() != coveredPoints(c, numPoints) * c.dim)
      return "a channel's values must hold dim floats per covered point";
    if (c.encoding != SpzChannelEncoding::Float32 &&
        !std::all_of(c.values.begin(), c.values.end(), [](float v) { return std::isfinite(v); })) {
      return "quantized channel values must be finite";
    }
    if (isStep(c)) {
      for (size_t i = 0; i < c.values.size(); ++i) {
        if (std::abs(static_cast<double>(c.values[i]) / c.steps[i % c.dim]) >
            static_cast<double>(std::numeric_limits<int32_t>::max())) {
          return "a step channel value is more than 2^31 steps from 0";
        }
      }
    }
    if (c.target == SpzTemporalTarget::Opacity && isSupported(c) &&
        c.dim == kWindowComponentsWithInterval) {
      for (size_t row = 0; row < c.values.size(); row += c.dim) {
        if (!std::isfinite(c.values[row + kWindowLo]) || !std::isfinite(c.values[row + kWindowHi]))
          return "window intervals must be finite";
      }
    }
    if (c.target != SpzTemporalTarget::Transform || !isSupported(c)) continue;
    const size_t influences = c.dim / kSkinValuesPerInfluence;
    for (size_t row = 0; row < c.values.size(); row += c.dim) {
      for (size_t j = 0; j < influences; ++j) {
        const float index = std::round(c.values[row + j]);
        const float weight = c.values[row + influences + j];
        if (!(index >= 0.0f && index < static_cast<float>(r.driverCount)))
          return "a skin index is outside the record's drivers";
        if (!std::isfinite(weight) || weight < 0.0f) return "skin weights must be finite and >= 0";
      }
    }
  }
  return {};
}

// Returns why the record can't be saved with, or applied to, a cloud of numPoints points, or an
// empty string if it can.
std::string recordProblem(const SpzExtensionTemporalNiantic& r, int32_t numPoints) {
  std::string problem = fieldsProblem(r);
  if (problem.empty()) problem = driverTransformsProblem(r);
  if (problem.empty()) problem = channelDataProblem(r, numPoints);
  return problem;
}

// ---- Record encoding ----

// Per-component bounds of a Range channel, computed when saving and read when loading.
struct ChannelBounds {
  std::vector<float> min;
  std::vector<float> max;
};

// Encodes a channel's values for points covered points, or returns std::nullopt (after logging
// why). Preconditions: the record passed recordProblem; bounds are the channel's Range bounds.
std::optional<std::vector<uint8_t>> encodeChannel(const SpzTemporalChannel& c, size_t points,
                                                  const ChannelBounds& bounds) {
  if (isStep(c)) {
    std::vector<uint8_t> out;
    if (!encodeComponentMajorSteps(c.values, points, c.dim, c.steps, &out)) {
      SpzLog("[SPZ ERROR] Temporal record: a step channel value can't be encoded");
      return std::nullopt;
    }
    return out;
  }
  return encodeComponentMajor(c.values, points, c.dim, isRange(c) ? c.bits : kFloat32Bits,
                              bounds.min, bounds.max);
}

// Decodes a channel's values. Returns false if a Range level exceeds its bits. Preconditions: the
// encoding is known; data holds channelBytes(c, points) bytes.
bool decodeChannel(const uint8_t* data, size_t points, const ChannelBounds& bounds,
                   SpzTemporalChannel* c) {
  if (isStep(*c)) {
    decodeComponentMajorSteps(data, points, c->dim, c->steps, &c->values);
    return true;
  }
  return decodeComponentMajor(data, points, c->dim, isRange(*c) ? c->bits : kFloat32Bits,
                              bounds.min, bounds.max, &c->values);
}

// Driver data in the file is component-major per driver ((driver * 7 + component) * K + keyframe),
// so each component's values over time sit together and compress well.
void appendDrivers(const SpzExtensionTemporalNiantic& r, std::vector<uint8_t>* out) {
  const size_t keyframes = r.driverKeyframeCount();
  const size_t start = out->size();
  out->resize(start + r.driverTransforms.size() * sizeof(float));
  for (size_t d = 0; d < r.driverCount; ++d) {
    for (size_t c = 0; c < kDriverTransformFloats; ++c) {
      for (size_t k = 0; k < keyframes; ++k) {
        const float v = r.driverTransforms[(d * keyframes + k) * kDriverTransformFloats + c];
        std::memcpy(&(*out)[start + ((d * kDriverTransformFloats + c) * keyframes + k) * sizeof(float)],
                    &v, sizeof(float));
      }
    }
  }
}

void decodeDrivers(const uint8_t* data, SpzExtensionTemporalNiantic* r) {
  const size_t keyframes = r->driverKeyframeCount();
  r->driverTransforms.assign(size_t{r->driverCount} * keyframes * kDriverTransformFloats, 0.0f);
  for (size_t d = 0; d < r->driverCount; ++d) {
    for (size_t c = 0; c < kDriverTransformFloats; ++c) {
      for (size_t k = 0; k < keyframes; ++k) {
        std::memcpy(&r->driverTransforms[(d * keyframes + k) * kDriverTransformFloats + c],
                    data + ((d * kDriverTransformFloats + c) * keyframes + k) * sizeof(float),
                    sizeof(float));
      }
    }
  }
}

uint64_t driverBytes(const SpzExtensionTemporalNiantic& r) {
  return uint64_t{r.driverCount} * r.driverKeyframeCount() * kDriverTransformFloats * sizeof(float);
}

// The channel descriptor fields that follow the target, type, flags, dim and point range.
void writeEncoding(std::ostream& os, const SpzTemporalChannel& c, const ChannelBounds& bounds) {
  writeValue(os, static_cast<uint8_t>(c.encoding));
  if (isRange(c)) {
    writeValue(os, c.bits);
    os.write(reinterpret_cast<const char*>(bounds.min.data()), c.dim * sizeof(float));
    os.write(reinterpret_cast<const char*>(bounds.max.data()), c.dim * sizeof(float));
  } else if (isStep(c)) {
    os.write(reinterpret_cast<const char*>(c.steps.data()), c.dim * sizeof(float));
  }
}

void writeName(std::ostream& os, const std::string& name) {
  writeValue(os, static_cast<uint16_t>(name.size()));
  os.write(name.data(), static_cast<std::streamsize>(name.size()));
}

void writeRecord(std::ostream& os, const SpzExtensionTemporalNiantic& r,
                 const std::vector<ChannelBounds>& bounds, const std::vector<uint64_t>& channelSizes,
                 const SpzTemporalTimeIndex& timeIndex,
                 const std::vector<SpzFrameLocation>& frames) {
  SpzLog("[SPZ] Writing extension: TemporalNiantic (clip %.3f-%.3fs, %zu channels, %u drivers)",
         r.clipStart, r.clipEnd, r.channels.size(), static_cast<unsigned>(r.driverCount));
  std::ostringstream payload;
  writeValue(payload, r.clipStart);
  writeValue(payload, r.clipEnd);
  writeValue(payload, r.posterTime);
  writeValue(payload, static_cast<uint8_t>(r.loop));
  writeValue(payload, uint8_t{0});  // reserved
  writeValue(payload, r.nominalFps);
  writeValue(payload, r.opacityCutoff);
  writeValue(payload, static_cast<uint8_t>(r.channels.size()));
  for (size_t i = 0; i < r.channels.size(); ++i) {
    const SpzTemporalChannel& c = r.channels[i];
    std::ostringstream descriptor;
    writeValue(descriptor, static_cast<uint8_t>(c.target));
    writeValue(descriptor, static_cast<uint8_t>(c.type));
    writeValue(descriptor, static_cast<uint8_t>(c.flags));
    writeValue(descriptor, c.dim);
    // pointEnd 0 means every point, since a real range never ends at 0.
    writeValue(descriptor, c.pointRange ? c.pointRange->begin : uint32_t{0});
    writeValue(descriptor, c.pointRange ? c.pointRange->end : uint32_t{0});
    // The data's size lets readers step over channels whose encoding they don't know.
    writeValue(descriptor, channelSizes[i]);
    writeEncoding(descriptor, c, bounds[i]);
    // Length-prefixed, so readers can skip descriptors with encoding parameters they don't know.
    const std::string bytes = descriptor.str();
    writeValue(payload, static_cast<uint32_t>(bytes.size()));
    payload.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  }
  writeValue(payload, r.driverCount);
  writeValue(payload, static_cast<uint32_t>(r.driverKeyframeCount()));
  for (double t : r.driverKeyframeTimes) writeValue(payload, t);
  writeValue(payload, static_cast<uint32_t>(r.driverNames.size()));
  for (const std::string& name : r.driverNames) writeName(payload, name);
  writeFrameTable(payload, frames);
  // The time index, for renderers; readers that predate it ignore bytes after the frame table.
  writeValue(payload, static_cast<uint32_t>(timeIndex.bins.size()));
  for (double t : timeIndex.boundaries) writeValue(payload, t);
  for (const SpzPointRange& bin : timeIndex.bins) {
    writeValue(payload, bin.begin);
    writeValue(payload, bin.end);
  }
  writeExtensionRecord(os, SpzExtensionType::SPZ_NIANTIC_temporal, payload.str());
}

// ---- Evaluation ----

struct RigidTransform {
  Vec3 translation;
  Quat rotation;  // unit
};

// Each driver's transform at clip time t: keyframes are interpolated linearly (translation) and by
// slerp (rotation), and held constant before the first and after the last.
std::vector<RigidTransform> sampleDrivers(const SpzExtensionTemporalNiantic& r, double t) {
  const std::vector<double>& times = r.driverKeyframeTimes;
  const size_t keyframes = times.size();
  size_t k0 = 0;
  size_t k1 = 0;
  float u = 0.0f;
  if (t >= times.back()) {
    k0 = k1 = keyframes - 1;
  } else if (t > times.front()) {
    k1 = static_cast<size_t>(std::upper_bound(times.begin(), times.end(), t) - times.begin());
    k0 = k1 - 1;
    u = static_cast<float>((t - times[k0]) / (times[k1] - times[k0]));
  }
  std::vector<RigidTransform> out(r.driverCount);
  for (size_t d = 0; d < r.driverCount; ++d) {
    const float* a = &r.driverTransforms[(d * keyframes + k0) * kDriverTransformFloats];
    const float* b = &r.driverTransforms[(d * keyframes + k1) * kDriverTransformFloats];
    out[d].translation = {a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1]),
                          a[2] + u * (b[2] - a[2])};
    out[d].rotation = slerp(normalize4(quatAt(a + kDriverRotation)),
                            normalize4(quatAt(b + kDriverRotation)), u);
  }
  return out;
}

// How many sigmas from its center a window stays at or above opacityCutoff; infinite for cutoff 0.
double windowReach(float opacityCutoff) {
  return opacityCutoff > 0.0f ? std::sqrt(-2.0 * std::log(static_cast<double>(opacityCutoff)))
                              : std::numeric_limits<double>::infinity();
}

// When a point is at least opacityCutoff opaque: from lo, until hi. The clip and intervals end
// before their end time ([lo, hi)); a window's reach includes its end (opacity equals the cutoff
// there).
struct VisibleSpan {
  double first;
  double second;
  bool includesEnd;
};

// The span during which a point with window row (nullptr if none) of windowDim components is at
// least opacityCutoff opaque, within the clip; std::nullopt if it never is. Shared by evaluation,
// the time index and the builder, so writers and readers agree.
std::optional<VisibleSpan> visibleSpanFor(const float* window, size_t windowDim,
                                          const SpzExtensionTemporalNiantic& clip) {
  VisibleSpan span{clip.clipStart, clip.clipEnd, false};
  const auto endAt = [&span](double hi, bool includesEnd) {
    if (hi < span.second || (hi == span.second && !includesEnd)) {
      span.second = hi;
      span.includesEnd = includesEnd;
    }
  };
  if (window) {
    const double reach = windowReach(clip.opacityCutoff) *
                         std::exp(static_cast<double>(window[kWindowLogSigma]));
    if (std::isfinite(reach)) {
      span.first = std::max(span.first, window[kWindowCenter] - reach);
      endAt(window[kWindowCenter] + reach, true);
    }
    if (windowDim == kWindowComponentsWithInterval) {
      span.first = std::max(span.first, static_cast<double>(window[kWindowLo]));
      endAt(window[kWindowHi], false);
    }
  }
  if (!(span.second > span.first)) return std::nullopt;
  return span;
}

// The time frame for a point with window row (nullptr if none) of windowDim components.
SpzTemporalTimeFrame timeFrameFor(const float* window, size_t windowDim,
                                  const SpzExtensionTemporalNiantic& clip) {
  const auto span = visibleSpanFor(window, windowDim, clip);
  if (!span) return {clip.posterTime, clip.clipEnd - clip.clipStart};
  return {std::clamp(clip.posterTime, span->first, span->second), span->second - span->first};
}

// The row of channel c holding point's values, or nullptr if c doesn't cover point.
const float* channelRow(const SpzTemporalChannel& c, size_t point) {
  const size_t first = c.pointRange ? c.pointRange->begin : 0;
  if (point < first) return nullptr;
  const size_t row = point - first;
  if ((row + 1) * c.dim > c.values.size()) return nullptr;
  return &c.values[row * c.dim];
}

const float* windowRow(const SpzTemporalModel& model, size_t point) {
  return model.opacity ? channelRow(*model.opacity, point) : nullptr;
}

size_t windowDim(const SpzTemporalModel& model) { return model.opacity ? model.opacity->dim : 0; }

using TimePowers = std::array<float, kMaxPolynomialDegree + 1>;

// u^0..u^3 for the point's normalized time u = (t - origin) / span.
TimePowers timePowers(const SpzTemporalModel& model, size_t point, double t) {
  const SpzTemporalTimeFrame frame = timeFrameFor(windowRow(model, point), windowDim(model),
                                                  *model.clip);
  const float u = static_cast<float>((t - frame.origin) / frame.span);
  return {1.0f, u, u * u, u * u * u};
}

// The points an evaluation computes; its output arrays hold them from index 0.
struct EvalSpan {
  size_t begin;
  size_t end;
};

// The points both channel c and span cover.
std::pair<size_t, size_t> overlap(const SpzTemporalChannel& c, int32_t numPoints, EvalSpan span) {
  const auto [first, last] = coveredRange(c, numPoints);
  return {std::max(first, span.begin), std::min(last, span.end)};
}

// Adds a polynomial channel's terms to values (components floats per point of span).
void addPolynomial(const SpzTemporalChannel& c, size_t components, const SpzTemporalModel& model,
                   double t, int32_t numPoints, EvalSpan span, std::vector<float>* values) {
  const size_t first = c.pointRange ? c.pointRange->begin : 0;
  const auto [from, to] = overlap(c, numPoints, span);
  const size_t degree = c.dim / components;
  for (size_t p = from; p < to; ++p) {
    const TimePowers powers = timePowers(model, p, t);
    const float* terms = &c.values[(p - first) * c.dim];
    float* out = &(*values)[(p - span.begin) * components];
    for (size_t k = 1; k <= degree; ++k) {
      for (size_t i = 0; i < components; ++i) out[i] += terms[(k - 1) * components + i] * powers[k];
    }
  }
}

void applyRotationPolynomial(const SpzTemporalChannel& c, const SpzTemporalModel& model, double t,
                             int32_t numPoints, EvalSpan span, SpzTemporalFrame* out) {
  addPolynomial(c, kQuaternionComponents, model, t, numPoints, span, &out->rotations);
  const auto [from, to] = overlap(c, numPoints, span);
  for (size_t p = from; p < to; ++p) {
    float* q = &out->rotations[(p - span.begin) * kQuaternionComponents];
    const Quat n = normalize4(quatAt(q));
    std::copy(n.begin(), n.end(), q);
  }
}

void applySkin(const SpzTemporalChannel& c, const std::vector<RigidTransform>& drivers,
               int32_t numPoints, EvalSpan span, SpzTemporalFrame* out) {
  const size_t first = c.pointRange ? c.pointRange->begin : 0;
  const auto [from, to] = overlap(c, numPoints, span);
  const size_t influences = c.dim / kSkinValuesPerInfluence;
  for (size_t p = from; p < to; ++p) {
    const float* row = &c.values[(p - first) * c.dim];
    float weightSum = 0.0f;
    for (size_t j = 0; j < influences; ++j) weightSum += row[influences + j];
    if (!(weightSum > 0.0f)) continue;  // unskinned: keeps its local state

    float* position = &out->positions[(p - span.begin) * kPositionComponents];
    const Vec3 local = {position[0], position[1], position[2]};
    Vec3 blendedPosition = {0.0f, 0.0f, 0.0f};
    Quat blendedRotation = {0.0f, 0.0f, 0.0f, 0.0f};
    const Quat* reference = nullptr;
    for (size_t j = 0; j < influences; ++j) {
      const float w = row[influences + j] / weightSum;
      if (w == 0.0f) continue;
      const RigidTransform& driver = drivers[static_cast<size_t>(std::round(row[j]))];
      const Vec3 moved = rotate(driver.rotation, local);
      for (size_t i = 0; i < kPositionComponents; ++i) {
        blendedPosition[i] += w * (moved[i] + driver.translation[i]);
      }
      // q and -q are the same rotation; align signs so they don't cancel in the blend.
      if (!reference) reference = &driver.rotation;
      const float sign = dot4(driver.rotation, *reference) < 0.0f ? -1.0f : 1.0f;
      for (size_t i = 0; i < kQuaternionComponents; ++i) {
        blendedRotation[i] += sign * w * driver.rotation[i];
      }
    }
    std::copy(blendedPosition.begin(), blendedPosition.end(), position);
    float* q = &out->rotations[(p - span.begin) * kQuaternionComponents];
    const Quat rotated = multiply(normalize4(blendedRotation), quatAt(q));
    std::copy(rotated.begin(), rotated.end(), q);
  }
}

// The opacity of a window row of windowDim components at clip time t.
double windowOpacity(const float* v, size_t windowDim, double t) {
  if (windowDim == kWindowComponentsWithInterval && !(t >= v[kWindowLo] && t < v[kWindowHi])) {
    return 0.0;
  }
  const double sigma = std::exp(static_cast<double>(v[kWindowLogSigma]));
  const double offset = (t - v[kWindowCenter]) / sigma;
  const double peak = std::clamp(static_cast<double>(v[kWindowPeak]), 0.0, 1.0);
  return peak * std::exp(-0.5 * offset * offset);
}

void applyWindow(const SpzTemporalChannel& c, double t, int32_t numPoints, EvalSpan span,
                 SpzTemporalFrame* out) {
  const size_t first = c.pointRange ? c.pointRange->begin : 0;
  const auto [from, to] = overlap(c, numPoints, span);
  for (size_t p = from; p < to; ++p) {
    out->alphas[p - span.begin] = alphaLogit(windowOpacity(&c.values[(p - first) * c.dim], c.dim, t));
  }
}

// Evaluates the points of span at clip time t into *out (sized for span).
void evaluateSpan(const GaussianCloud& cloud, const SpzTemporalModel& model, double t, EvalSpan span,
                  SpzTemporalFrame* out) {
  const auto slice = [&](const std::vector<float>& from, size_t components, std::vector<float>* to) {
    to->assign(from.begin() + static_cast<std::ptrdiff_t>(span.begin * components),
               from.begin() + static_cast<std::ptrdiff_t>(span.end * components));
  };
  slice(cloud.positions, kPositionComponents, &out->positions);
  slice(cloud.rotations, kQuaternionComponents, &out->rotations);
  slice(cloud.alphas, 1, &out->alphas);
  if (t < model.clip->clipStart || t >= model.clip->clipEnd) {
    std::fill(out->alphas.begin(), out->alphas.end(), alphaLogit(0.0));
    return;
  }
  const int32_t n = cloud.numPoints;
  if (model.position) {
    addPolynomial(*model.position, kPositionComponents, model, t, n, span, &out->positions);
  }
  if (model.rotation) applyRotationPolynomial(*model.rotation, model, t, n, span, out);
  if (model.skin) applySkin(*model.skin, sampleDrivers(*model.clip, t), n, span, out);
  if (model.opacity) applyWindow(*model.opacity, t, n, span, out);
}

// A non-owning model of record, for computing its time index while saving it.
SpzTemporalModel unownedModel(const SpzExtensionTemporalNiantic& record) {
  SpzTemporalModel model;
  // Aliasing constructor with an empty owner: points at record without owning it.
  model.clip = std::shared_ptr<const SpzExtensionTemporalNiantic>(std::shared_ptr<void>(), &record);
  for (const SpzTemporalChannel& c : record.channels) {
    if (c.target == SpzTemporalTarget::Opacity && isSupported(c)) model.opacity = &c;
  }
  return model;
}

// ---- Builder ----

// C(j, k) for j, k <= kMaxPolynomialDegree.
constexpr std::array<std::array<double, 4>, 4> kBinomial = {{
  {1, 0, 0, 0},
  {1, 1, 0, 0},
  {1, 2, 1, 0},
  {1, 3, 3, 1},
}};

// Rewrites the polynomial value + Σ_{j=1..degree} coeffs[j-1] (t - c)^j about the new origin
// c + delta: value becomes the polynomial at c + delta, and coeffs the coefficients in powers of
// (t - c - delta). components is the number of floats per term.
void recenter(float* value, float* coeffs, size_t degree, size_t components, double delta) {
  const std::array<double, kMaxPolynomialDegree + 1> powers = {1.0, delta, delta * delta,
                                                               delta * delta * delta};
  for (size_t c = 0; c < components; ++c) {
    std::array<double, kMaxPolynomialDegree + 1> a = {value[c], 0, 0, 0};
    for (size_t j = 1; j <= degree; ++j) a[j] = coeffs[(j - 1) * components + c];
    // a'_k = Σ_{j>=k} C(j, k) a_j delta^(j-k)
    for (size_t k = 0; k <= degree; ++k) {
      double sum = 0.0;
      for (size_t j = k; j <= degree; ++j) sum += kBinomial[j][k] * a[j] * powers[j - k];
      if (k == 0) {
        value[c] = static_cast<float>(sum);
      } else {
        coeffs[(k - 1) * components + c] = static_cast<float>(sum);
      }
    }
  }
}

bool allFinite(const std::vector<float>& v) {
  return std::all_of(v.begin(), v.end(), [](float x) { return std::isfinite(x); });
}

// Returns why source can't be built, or an empty string if it can.
std::string sourceProblem(const SpzTemporalSource& s) {
  const GaussianCloud& g = s.cloud;
  const size_t n = static_cast<size_t>(std::max(g.numPoints, 0));
  if (n == 0) return "the cloud has no points";
  const size_t shFloats = static_cast<size_t>(dimForDegree(g.shDegree)) * 3;
  if (g.positions.size() != n * 3 || g.scales.size() != n * 3 || g.rotations.size() != n * 4 ||
      g.alphas.size() != n || g.colors.size() != n * 3 || g.sh.size() != n * shFloats) {
    return "the cloud's attribute sizes don't match numPoints";
  }
  for (const auto& ext : g.extensions) {
    if (ext->extensionType == SpzExtensionType::SPZ_NIANTIC_custom_attribute ||
        ext->extensionType == SpzExtensionType::SPZ_NIANTIC_temporal) {
      return "the cloud must not carry custom attributes or a temporal record";
    }
  }
  if (!s.tCenter.empty() && (s.tCenter.size() != n || !allFinite(s.tCenter)))
    return "tCenter must be empty or hold numPoints finite times";
  if (!s.tLogSigma.empty() && (s.tLogSigma.size() != n || !allFinite(s.tLogSigma)))
    return "tLogSigma must be empty or hold numPoints finite values";
  if (!s.tInterval.empty() && (s.tInterval.size() != n * 2 || !allFinite(s.tInterval)))
    return "tInterval must be empty or hold numPoints * 2 finite times";
  if (s.positionDegree > kMaxPolynomialDegree || s.rotationDegree > kMaxPolynomialDegree)
    return "polynomial degrees must be 0-3";
  if (s.positionCoefficients.size() != n * kPositionComponents * s.positionDegree ||
      !allFinite(s.positionCoefficients)) {
    return "positionCoefficients must hold numPoints * 3 * positionDegree finite floats";
  }
  if (s.rotationCoefficients.size() != n * kQuaternionComponents * s.rotationDegree ||
      !allFinite(s.rotationCoefficients)) {
    return "rotationCoefficients must hold numPoints * 4 * rotationDegree finite floats";
  }
  const SpzTemporalPrecision& e = s.precision;
  if (!validSteps({e.position, e.rotation, e.time, e.logSigma, e.opacity}))
    return "precision bounds must be finite and positive";
  SpzExtensionTemporalNiantic clip;
  clip.clipStart = s.clipStart;
  clip.clipEnd = s.clipEnd;
  clip.posterTime = s.posterTime;
  clip.loop = s.loop;
  clip.nominalFps = s.nominalFps;
  clip.opacityCutoff = s.opacityCutoff;
  return fieldsProblem(clip);
}

// Copies point src of from into point dst of to, for every core attribute.
void copyPoint(const GaussianCloud& from, size_t src, size_t dst, size_t shFloats,
               GaussianCloud* to) {
  std::copy_n(&from.positions[src * 3], 3, &to->positions[dst * 3]);
  std::copy_n(&from.scales[src * 3], 3, &to->scales[dst * 3]);
  std::copy_n(&from.rotations[src * 4], 4, &to->rotations[dst * 4]);
  to->alphas[dst] = from.alphas[src];
  std::copy_n(&from.colors[src * 3], 3, &to->colors[dst * 3]);
  std::copy_n(from.sh.begin() + static_cast<std::ptrdiff_t>(src * shFloats), shFloats,
              to->sh.begin() + static_cast<std::ptrdiff_t>(dst * shFloats));
}
}  // namespace

// ---- SpzExtensionTemporalNiantic ----

SpzExtensionTemporalNiantic::SpzExtensionTemporalNiantic()
    : SpzExtensionBase(SpzExtensionType::SPZ_NIANTIC_temporal) {}

SpzExtensionType SpzExtensionTemporalNiantic::type() {
  return SpzExtensionType::SPZ_NIANTIC_temporal;
}

SpzExtensionBase* SpzExtensionTemporalNiantic::copyAsRawData() const {
  return new SpzExtensionTemporalNiantic(*this);
}

uint32_t SpzExtensionTemporalNiantic::payloadBytes() const { return 0; }

void SpzExtensionTemporalNiantic::write(std::ostream& /*os*/) const {
  SpzLog("[SPZ ERROR] Temporal record not written: it can only be saved by saveSpz in SPZ "
         "version 4 files");
}

const SpzTemporalChannel* SpzExtensionTemporalNiantic::findChannel(SpzTemporalTarget target) const {
  for (const SpzTemporalChannel& c : channels) {
    if (c.target == target) return &c;
  }
  return nullptr;
}

std::optional<SpzFrameWritePlan> SpzExtensionTemporalNiantic::planFrameWrites(
    int32_t numPoints) const {
  const std::string problem = recordProblem(*this, numPoints);
  if (!problem.empty()) {
    SpzLog("[SPZ ERROR] Temporal record: %s", problem.c_str());
    return std::nullopt;
  }
  // Data: the drivers, then each channel's values, in channel order.
  std::vector<uint8_t> data;
  appendDrivers(*this, &data);
  std::vector<ChannelBounds> bounds(channels.size());
  std::vector<uint64_t> channelSizes;
  for (size_t i = 0; i < channels.size(); ++i) {
    const SpzTemporalChannel& c = channels[i];
    const size_t points = coveredPoints(c, numPoints);
    if (!channelBytes(c, points)) {
      SpzLog("[SPZ ERROR] Temporal record: a channel is too large for this platform");
      return std::nullopt;
    }
    if (isRange(c)) componentRanges(c.values, c.dim, &bounds[i].min, &bounds[i].max);
    const std::optional<std::vector<uint8_t>> encoded = encodeChannel(c, points, bounds[i]);
    if (!encoded) return std::nullopt;
    channelSizes.push_back(encoded->size());
    data.insert(data.end(), encoded->begin(), encoded->end());
  }
  // Index the values readers will decode, so the index holds for them after quantization.
  SpzExtensionTemporalNiantic decoded(*this);
  size_t offset = static_cast<size_t>(driverBytes(*this));
  for (size_t i = 0; i < decoded.channels.size(); ++i) {
    SpzTemporalChannel& c = decoded.channels[i];
    if (!decodeChannel(data.data() + offset, coveredPoints(c, numPoints), bounds[i], &c)) {
      SpzLog("[SPZ ERROR] Temporal record: a channel doesn't decode back");
      return std::nullopt;
    }
    offset += static_cast<size_t>(channelSizes[i]);
  }
  const SpzTemporalModel model = unownedModel(decoded);
  const SpzTemporalTimeIndex timeIndex =
    computeTemporalTimeIndex(model, numPoints, temporalTimeBinBoundaries(model));
  SpzExtensionTemporalNiantic copy(*this);
  return SpzFrameWritePlan{std::move(data),
                           [copy, bounds, channelSizes, timeIndex](
                             std::ostream& os, const std::vector<SpzFrameLocation>& frames) {
                             writeRecord(os, copy, bounds, channelSizes, timeIndex, frames);
                           }};
}

std::string SpzExtensionTemporalNiantic::uniquenessKey() const { return kRecordKey; }

SpzExtensionBasePtr SpzExtensionTemporalNiantic::withConvertedCoordinates(
    const CoordinateConverter& c) const {
  const bool positionsChange = !isIdentityPositionTransform(c);
  const bool rotationsChange = !isIdentityRotationTransform(c);
  if (!positionsChange && !rotationsChange) return nullptr;
  auto converted = std::make_shared<SpzExtensionTemporalNiantic>(*this);
  for (SpzTemporalChannel& channel : converted->channels) {
    // Polynomial terms transform like the core values they're added to.
    if (isPositionPolynomial(channel) && positionsChange) {
      for (size_t i = 0; i + kPositionComponents <= channel.values.size();
           i += kPositionComponents) {
        convertDirection(c, &channel.values[i]);
      }
      unifyStepsPerVector(kPositionComponents, &channel);
    } else if (isRotationPolynomial(channel) && rotationsChange) {
      for (size_t i = 0; i + kQuaternionComponents <= channel.values.size();
           i += kQuaternionComponents) {
        convertQuaternion(c, &channel.values[i]);
      }
      unifyStepsPerVector(kQuaternionComponents, &channel);
    }
  }
  if (driverCount == 0 || !positionsChange) return converted;
  // The conversion as a matrix: its columns are the converted axes. It is orthogonal (a rotation
  // and/or axis flips), so a driver's rigid transform (R, T) becomes (C R Cᵀ, C T).
  Mat3 conversion{};
  for (size_t axis = 0; axis < kPositionComponents; ++axis) {
    float e[kPositionComponents] = {0.0f, 0.0f, 0.0f};
    e[axis] = 1.0f;
    convertDirection(c, e);
    for (size_t row = 0; row < kPositionComponents; ++row) conversion[row][axis] = e[row];
  }
  const Mat3 inverse = transposed(conversion);
  for (size_t i = 0; i < converted->driverTransforms.size(); i += kDriverTransformFloats) {
    float* f = &converted->driverTransforms[i];
    convertDirection(c, f);
    float* r = f + kDriverRotation;
    const Quat q =
      fromMatrix(multiply(multiply(conversion, toMatrix(normalize4(quatAt(r)))), inverse));
    std::copy(q.begin(), q.end(), r);
  }
  return converted;
}

SpzExtensionBasePtr SpzExtensionTemporalNiantic::withRotationFactors(
    const std::vector<float>& factors) const {
  const SpzTemporalChannel* rotation = findChannel(SpzTemporalTarget::Rotation);
  if (!rotation || !isRotationPolynomial(*rotation)) return nullptr;
  auto scaled = std::make_shared<SpzExtensionTemporalNiantic>(*this);
  for (SpzTemporalChannel& channel : scaled->channels) {
    if (!isRotationPolynomial(channel) || channel.dim == 0) continue;
    const size_t first = channel.pointRange ? channel.pointRange->begin : 0;
    const size_t points = channel.values.size() / channel.dim;
    for (size_t p = 0; p < points && first + p < factors.size(); ++p) {
      for (size_t i = 0; i < channel.dim; ++i) channel.values[p * channel.dim + i] *= factors[first + p];
    }
  }
  return scaled;
}

std::optional<std::shared_ptr<SpzExtensionBase>> SpzExtensionTemporalNiantic::tryReadFromPly(
    std::istream& /*in*/, const std::unordered_set<std::string>& /*elementNames*/) const {
  return std::nullopt;
}

void SpzExtensionTemporalNiantic::writePlyHeader(std::ostream& /*out*/) const {
  SpzLog("[SPZ WARNING] PLY can't store an animated clip: saved the poster frame only");
}

void SpzExtensionTemporalNiantic::writePlyData(std::ostream& /*out*/) const {}

namespace {
// A channel as its descriptor describes it, with what's needed to find and decode its data.
struct ParsedChannel {
  SpzTemporalChannel channel;
  ChannelBounds bounds;
  uint64_t dataBytes = 0;
  // Unset for an optional channel whose encoding this reader doesn't know: its data is skipped.
  bool known = true;
};

// Reads one length-prefixed channel descriptor, or returns std::nullopt (after logging why) if it
// is malformed or a required channel uses an unknown encoding.
std::optional<ParsedChannel> readChannel(std::istream& is) {
  using detail::readExact;
  uint32_t descriptorBytes = 0;
  if (!readExact(is, descriptorBytes)) {
    logStatic("truncated channel");
    return std::nullopt;
  }
  const std::streamoff remaining = detail::remainingBytes(is);
  if (remaining < 0 || descriptorBytes > static_cast<uint64_t>(remaining)) {
    logStatic("truncated channel");
    return std::nullopt;
  }
  std::string bytes(descriptorBytes, '\0');
  is.read(bytes.data(), descriptorBytes);
  std::istringstream d(bytes);

  ParsedChannel parsed;
  SpzTemporalChannel& c = parsed.channel;
  uint32_t begin = 0;
  uint32_t end = 0;
  uint8_t encoding = 0;
  if (!readExact(d, c.target) || !readExact(d, c.type) || !readExact(d, c.flags) ||
      !readExact(d, c.dim) || !readExact(d, begin) || !readExact(d, end) ||
      !readExact(d, parsed.dataBytes) || !readExact(d, encoding)) {
    logStatic("truncated channel");
    return std::nullopt;
  }
  if (end != 0) c.pointRange = SpzPointRange{begin, end};
  if (!isKnownEncoding(encoding)) {
    if (isRequired(c)) {
      logStatic("a required channel uses an unknown encoding");
      return std::nullopt;
    }
    SpzLog("[SPZ WARNING] Temporal record: skipped an optional channel with an unknown encoding");
    parsed.known = false;
    return parsed;
  }
  c.encoding = static_cast<SpzChannelEncoding>(encoding);
  if (isRange(c)) {
    parsed.bounds.min.resize(c.dim);
    parsed.bounds.max.resize(c.dim);
    if (!readExact(d, c.bits) ||
        !d.read(reinterpret_cast<char*>(parsed.bounds.min.data()), c.dim * sizeof(float)) ||
        !d.read(reinterpret_cast<char*>(parsed.bounds.max.data()), c.dim * sizeof(float))) {
      logStatic("truncated channel");
      return std::nullopt;
    }
    for (size_t k = 0; k < c.dim; ++k) {
      if (!std::isfinite(parsed.bounds.min[k]) || !std::isfinite(parsed.bounds.max[k]) ||
          parsed.bounds.min[k] > parsed.bounds.max[k]) {
        logStatic("invalid channel range");
        return std::nullopt;
      }
    }
  } else if (isStep(c)) {
    c.steps.resize(c.dim);
    if (!d.read(reinterpret_cast<char*>(c.steps.data()), c.dim * sizeof(float))) {
      logStatic("truncated channel");
      return std::nullopt;
    }
  }
  // Bytes after the encoding parameters are ignored: later versions may append fields.
  return parsed;
}

bool readName(std::istream& is, std::string* name) {
  uint16_t length = 0;
  if (!detail::readExact(is, length)) return false;
  name->resize(length);
  return static_cast<bool>(is.read(name->data(), length));
}
}  // namespace

std::optional<SpzFrameReadPlan> parseTemporalNiantic(std::istream& is) {
  using detail::readExact;
  auto r = std::make_shared<SpzExtensionTemporalNiantic>();
  uint8_t loop = 0;
  uint8_t reserved = 0;
  uint8_t channelCount = 0;
  if (!readExact(is, r->clipStart) || !readExact(is, r->clipEnd) ||
      !readExact(is, r->posterTime) || !readExact(is, loop) || !readExact(is, reserved) ||
      !readExact(is, r->nominalFps) || !readExact(is, r->opacityCutoff) ||
      !readExact(is, channelCount)) {
    logStatic("truncated record");
    return std::nullopt;
  }
  // The loop mode is only a playback hint, so an unknown one falls back to playing once.
  r->loop = loop <= kMaxLoop ? static_cast<SpzTemporalLoop>(loop) : SpzTemporalLoop::Once;
  std::vector<ParsedChannel> channels;
  for (uint8_t i = 0; i < channelCount; ++i) {
    std::optional<ParsedChannel> channel = readChannel(is);
    if (!channel) return std::nullopt;
    channels.push_back(std::move(*channel));
  }
  for (const ParsedChannel& c : channels) {
    if (c.known) r->channels.push_back(c.channel);
  }
  uint32_t keyframes = 0;
  if (!readExact(is, r->driverCount) || !readExact(is, keyframes)) {
    logStatic("truncated drivers");
    return std::nullopt;
  }
  // Checked before allocating, so a count that lies can't drive a large allocation.
  const std::streamoff remaining = detail::remainingBytes(is);
  if (remaining < 0 || uint64_t{keyframes} * sizeof(double) > static_cast<uint64_t>(remaining)) {
    logStatic("truncated driver keyframe times");
    return std::nullopt;
  }
  r->driverKeyframeTimes.resize(keyframes);
  for (double& t : r->driverKeyframeTimes) readExact(is, t);
  uint32_t nameCount = 0;
  if (!readExact(is, nameCount) || nameCount > r->driverCount) {
    logStatic("invalid driver names");
    return std::nullopt;
  }
  r->driverNames.resize(nameCount);
  for (std::string& name : r->driverNames) {
    if (!readName(is, &name)) {
      logStatic("truncated driver names");
      return std::nullopt;
    }
  }
  std::optional<std::vector<SpzFrameLocation>> frames = readFrameTable(is);
  if (!frames) {
    logStatic("truncated frame table");
    return std::nullopt;
  }
  // The time index follows; this library recomputes it from the data rather than trust it, and
  // ignores it and any later fields.
  const std::string problem = fieldsProblem(*r);
  if (!problem.empty()) {
    logStatic(problem);
    return std::nullopt;
  }
  SpzLog("[SPZ] Found extension: TemporalNiantic (clip %.3f-%.3fs, %zu channels, %u drivers)",
         r->clipStart, r->clipEnd, r->channels.size(), static_cast<unsigned>(r->driverCount));

  // Whether every known channel's declared data size matches its descriptor for numPoints points
  // (logged otherwise).
  const auto sizesMatch = [channels](int32_t numPoints) {
    for (const ParsedChannel& c : channels) {
      if (!c.known) continue;
      if (const char* rangeIssue = rangeProblem(c.channel, numPoints)) {
        logStatic(rangeIssue);
        return false;
      }
      const std::optional<size_t> bytes = channelBytes(c.channel, coveredPoints(c.channel, numPoints));
      if (!bytes || *bytes != c.dataBytes) {
        logStatic("a channel's data size doesn't match its descriptor");
        return false;
      }
    }
    return true;
  };
  return SpzFrameReadPlan{
    SpzExtensionType::SPZ_NIANTIC_temporal,
    std::move(*frames),
    std::string(kRecordKey),
    [r, channels, sizesMatch](int32_t numPoints) -> std::optional<uint64_t> {
      if (!sizesMatch(numPoints)) return std::nullopt;
      uint64_t total = driverBytes(*r);
      for (const ParsedChannel& c : channels) {
        if (c.dataBytes > std::numeric_limits<uint64_t>::max() - total) return std::nullopt;
        total += c.dataBytes;
      }
      return total;
    },
    [r, channels](int32_t numPoints, std::vector<uint8_t>&& data) -> SpzExtensionBasePtr {
      auto resolved = std::make_shared<SpzExtensionTemporalNiantic>(*r);
      resolved->channels.clear();
      decodeDrivers(data.data(), resolved.get());
      size_t offset = static_cast<size_t>(driverBytes(*r));
      for (const ParsedChannel& parsed : channels) {
        if (parsed.known) {
          SpzTemporalChannel c = parsed.channel;
          if (!decodeChannel(data.data() + offset, coveredPoints(c, numPoints), parsed.bounds, &c)) {
            logStatic("a channel's quantized level is out of range");
            return nullptr;
          }
          resolved->channels.push_back(std::move(c));
        }
        offset += static_cast<size_t>(parsed.dataBytes);
      }
      const std::string resolvedProblem = recordProblem(*resolved, numPoints);
      if (!resolvedProblem.empty()) {
        logStatic(resolvedProblem);
        return nullptr;
      }
      return resolved;
    }};
}

// ---- Binding and evaluation ----

std::optional<SpzTemporalModel> bindTemporalModel(const GaussianCloud& cloud) {
  auto record = findExtensionByType<SpzExtensionTemporalNiantic>(cloud.extensions);
  if (!record) return std::nullopt;
  const std::string problem = recordProblem(*record, cloud.numPoints);
  if (!problem.empty()) {
    logStatic(problem);
    return std::nullopt;
  }
  SpzTemporalModel model;
  model.clip = record;
  for (const SpzTemporalChannel& c : record->channels) {
    if (!isSupported(c)) continue;  // optional, or recordProblem would have failed
    switch (c.target) {
      case SpzTemporalTarget::Position: model.position = &c; break;
      case SpzTemporalTarget::Rotation: model.rotation = &c; break;
      case SpzTemporalTarget::Opacity: model.opacity = &c; break;
      case SpzTemporalTarget::Transform: model.skin = &c; break;
    }
  }
  return model;
}

void evaluateTemporalInto(const GaussianCloud& cloud, const SpzTemporalModel& model, double t,
                          SpzTemporalFrame* frame) {
  evaluateSpan(cloud, model, t, {0, static_cast<size_t>(cloud.numPoints)}, frame);
}

std::vector<double> uniformTimeBinBoundaries(double start, double end, size_t count) {
  std::vector<double> boundaries(count + 1);
  for (size_t b = 0; b <= count; ++b) {
    boundaries[b] = start + (end - start) * static_cast<double>(b) / static_cast<double>(count);
  }
  boundaries.back() = end;
  return boundaries;
}

std::vector<double> temporalTimeBinBoundaries(const SpzTemporalModel& model) {
  const SpzExtensionTemporalNiantic& clip = *model.clip;
  if (model.opacity && model.opacity->dim == kWindowComponentsWithInterval) {
    std::vector<double> bounds = {clip.clipStart, clip.clipEnd};
    const std::vector<float>& v = model.opacity->values;
    for (size_t row = 0; row + kWindowComponentsWithInterval <= v.size();
         row += kWindowComponentsWithInterval) {
      for (size_t k : {kWindowLo, kWindowHi}) {
        const double t = v[row + k];
        if (t > clip.clipStart && t < clip.clipEnd) bounds.push_back(t);
      }
    }
    std::sort(bounds.begin(), bounds.end());
    bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
    if (bounds.size() - 1 <= kMaxTimeBins) return bounds;
  }
  size_t count = kDefaultTimeBins;
  if (clip.nominalFps > 0.0f) {
    const double frames = std::ceil((clip.clipEnd - clip.clipStart) * clip.nominalFps);
    count = static_cast<size_t>(std::clamp(frames, 1.0, static_cast<double>(kMaxTimeBins)));
  }
  return uniformTimeBinBoundaries(clip.clipStart, clip.clipEnd, count);
}

SpzTemporalTimeIndex computeTemporalTimeIndex(const SpzTemporalModel& model, int32_t numPoints,
                                              const std::vector<double>& boundaries) {
  const SpzExtensionTemporalNiantic& clip = *model.clip;
  const size_t binCount = boundaries.size() - 1;
  SpzTemporalTimeIndex index;
  index.boundaries = boundaries;
  constexpr uint32_t kNone = std::numeric_limits<uint32_t>::max();
  std::vector<SpzPointRange> bins(binCount, SpzPointRange{kNone, 0});
  // Points visible throughout the clip (or never measured against a window) widen every bin.
  SpzPointRange everywhere{kNone, 0};
  const auto widen = [](SpzPointRange& range, uint32_t point) {
    range.begin = std::min(range.begin, point);
    range.end = std::max(range.end, point + 1);
  };
  const auto clampBin = [&](std::ptrdiff_t b) {
    return static_cast<size_t>(std::clamp<std::ptrdiff_t>(b, 0, static_cast<std::ptrdiff_t>(binCount) - 1));
  };
  // The bin holding t, and the last bin a span reaches: one ending exactly on a boundary,
  // exclusively, stops before it.
  const auto binOf = [&](double t) {
    return clampBin(std::upper_bound(boundaries.begin(), boundaries.end(), t) - boundaries.begin() - 1);
  };
  const auto lastBinOf = [&](const VisibleSpan& span) {
    const auto end = span.includesEnd
                       ? std::upper_bound(boundaries.begin(), boundaries.end(), span.second)
                       : std::lower_bound(boundaries.begin(), boundaries.end(), span.second);
    return clampBin(end - boundaries.begin() - 1);
  };
  for (uint32_t p = 0; p < static_cast<uint32_t>(numPoints); ++p) {
    const float* row = windowRow(model, p);
    if (!row) {
      widen(everywhere, p);
      continue;
    }
    const auto span = visibleSpanFor(row, windowDim(model), clip);
    if (!span) continue;
    const size_t first = binOf(span->first);
    const size_t last = lastBinOf(*span);
    if (first == 0 && last == binCount - 1) {
      widen(everywhere, p);
      continue;
    }
    for (size_t b = first; b <= last; ++b) widen(bins[b], p);
  }
  for (SpzPointRange& bin : bins) {
    if (everywhere.begin != kNone) {
      bin.begin = std::min(bin.begin, everywhere.begin);
      bin.end = std::max(bin.end, everywhere.end);
    }
    if (bin.begin == kNone) bin = {0, 0};
    index.maxBinPoints = std::max(index.maxBinPoints, bin.end - bin.begin);
  }
  index.bins = std::move(bins);
  return index;
}

SpzPointRange temporalCandidates(const SpzTemporalTimeIndex& index, double t) {
  const std::vector<double>& b = index.boundaries;
  if (index.bins.empty() || !(t >= b.front() && t < b.back())) return {0, 0};
  const size_t bin = static_cast<size_t>(std::upper_bound(b.begin(), b.end(), t) - b.begin() - 1);
  return index.bins[std::min(bin, index.bins.size() - 1)];
}

void evaluateVisibleTemporal(const GaussianCloud& cloud, const SpzTemporalModel& model,
                             const SpzTemporalTimeIndex& index, double t,
                             SpzTemporalVisibleFrame* frame) {
  const SpzPointRange candidates = temporalCandidates(index, t);
  evaluateSpan(cloud, model, t, {candidates.begin, candidates.end}, &frame->work);
  const float minAlpha = model.clip->opacityCutoff > 0.0f ? alphaLogit(model.clip->opacityCutoff)
                                                          : -std::numeric_limits<float>::infinity();
  frame->indices.clear();
  frame->positions.clear();
  frame->rotations.clear();
  frame->alphas.clear();
  const SpzTemporalFrame& work = frame->work;
  for (size_t i = 0; i < work.alphas.size(); ++i) {
    if (work.alphas[i] < minAlpha) continue;
    frame->indices.push_back(static_cast<uint32_t>(candidates.begin + i));
    frame->positions.insert(frame->positions.end(), &work.positions[i * kPositionComponents],
                            &work.positions[(i + 1) * kPositionComponents]);
    frame->rotations.insert(frame->rotations.end(), &work.rotations[i * kQuaternionComponents],
                            &work.rotations[(i + 1) * kQuaternionComponents]);
    frame->alphas.push_back(work.alphas[i]);
  }
}

SpzTemporalTimeFrame temporalTimeFrame(const SpzTemporalModel& model, size_t point) {
  return timeFrameFor(windowRow(model, point), model.opacity ? model.opacity->dim : 0, *model.clip);
}

GaussianCloud evaluateTemporal(const GaussianCloud& cloud, const SpzTemporalModel& model,
                               double t) {
  SpzTemporalFrame frame;
  evaluateTemporalInto(cloud, model, t, &frame);
  GaussianCloud out = cloud;
  out.positions = std::move(frame.positions);
  out.rotations = std::move(frame.rotations);
  out.alphas = std::move(frame.alphas);
  out.extensions.clear();
  for (const SpzExtensionBasePtr& ext : cloud.extensions) {
    if (ext != model.clip) out.extensions.push_back(ext);
  }
  return out;
}

double clipTime(const SpzExtensionTemporalNiantic& clip, double playbackSeconds) {
  const double duration = clip.clipEnd - clip.clipStart;
  const double lastInstant = std::nextafter(clip.clipEnd, clip.clipStart);
  const double s = std::max(playbackSeconds, 0.0);
  switch (clip.loop) {
    case SpzTemporalLoop::Loop:
      return std::min(clip.clipStart + std::fmod(s, duration), lastInstant);
    case SpzTemporalLoop::PingPong: {
      const double cycle = kPingPongPasses * duration;
      const double phase = std::fmod(s, cycle);
      const double offset = phase <= duration ? phase : cycle - phase;
      return std::min(clip.clipStart + offset, lastInstant);
    }
    case SpzTemporalLoop::Once:
    default:
      return std::min(clip.clipStart + s, lastInstant);
  }
}

// ---- Builder ----

std::optional<GaussianCloud> buildTemporalCloud(const SpzTemporalSource& source) {
  const std::string problem = sourceProblem(source);
  if (!problem.empty()) {
    SpzLog("[SPZ ERROR] buildTemporalCloud: %s", problem.c_str());
    return std::nullopt;
  }
  const GaussianCloud& g = source.cloud;
  const size_t n = static_cast<size_t>(g.numPoints);
  const size_t posDim = kPositionComponents * source.positionDegree;
  const size_t rotDim = kQuaternionComponents * source.rotationDegree;
  const bool hasInterval = !source.tInterval.empty();
  const bool windowed = !source.tLogSigma.empty() || hasInterval;
  const size_t windowDim = hasInterval ? kWindowComponentsWithInterval : kWindowComponents;
  const double duration = source.clipEnd - source.clipStart;
  const double alwaysVisibleSigma = kAlwaysVisibleSigmaRatio * duration;
  const float alwaysVisibleLogSigma = static_cast<float>(std::log(alwaysVisibleSigma));

  // A step of twice the bound keeps every quantized value within the bound.
  const SpzTemporalPrecision& precision = source.precision;
  const float positionStep = 2.0f * precision.position;
  const float rotationStep = 2.0f * precision.rotation;
  const float timeStep = 2.0f * precision.time;
  std::vector<float> windowSteps = {timeStep, 2.0f * precision.logSigma, 2.0f * precision.opacity,
                                    timeStep, timeStep};
  windowSteps.resize(windowDim);

  auto clip = std::make_shared<SpzExtensionTemporalNiantic>();
  clip->clipStart = source.clipStart;
  clip->clipEnd = source.clipEnd;
  clip->posterTime = source.posterTime;
  clip->loop = source.loop;
  clip->nominalFps = source.nominalFps;
  clip->opacityCutoff = source.opacityCutoff;

  // Work on copies in source order: each point's state at its origin goes into `poster`, and
  // coefficients are re-expressed in its time frame in place.
  GaussianCloud poster = g;
  std::vector<float> posCoeffs = source.positionCoefficients;
  std::vector<float> rotCoeffs = source.rotationCoefficients;
  std::vector<float> window(windowed ? n * windowDim : 0);
  std::vector<uint8_t> animated(n, 0);
  for (size_t i = 0; i < n; ++i) {
    const double center = source.tCenter.empty() ? source.posterTime : source.tCenter[i];
    float* row = windowed ? &window[i * windowDim] : nullptr;
    if (row) {
      row[kWindowCenter] = static_cast<float>(center);
      row[kWindowLogSigma] = source.tLogSigma.empty() ? alwaysVisibleLogSigma : source.tLogSigma[i];
      row[kWindowPeak] = sigmoid(g.alphas[i]);
      if (hasInterval) {
        row[kWindowLo] = source.tInterval[i * 2];
        row[kWindowHi] = source.tInterval[i * 2 + 1];
      }
      // Readers derive time frames from the stored values, so derive them from those too.
      for (size_t c = 0; c < windowDim; ++c) row[c] = quantizeToStep(row[c], windowSteps[c]);
    }
    const SpzTemporalTimeFrame frame = timeFrameFor(row, windowDim, *clip);
    const double delta = frame.origin - center;
    const std::array<double, kMaxPolynomialDegree + 1> spanPowers = {
      1.0, frame.span, frame.span * frame.span, frame.span * frame.span * frame.span};
    if (posDim > 0) {
      float* a = &posCoeffs[i * posDim];
      recenter(&poster.positions[i * 3], a, source.positionDegree, kPositionComponents, delta);
      for (size_t k = 1; k <= source.positionDegree; ++k)
        for (size_t c = 0; c < kPositionComponents; ++c)
          a[(k - 1) * kPositionComponents + c] *= static_cast<float>(spanPowers[k]);
    }
    if (rotDim > 0) {
      float* q = &poster.rotations[i * 4];
      float* b = &rotCoeffs[i * rotDim];
      recenter(q, b, source.rotationDegree, kQuaternionComponents, delta);
      // normalize(q + Σ b u^k) is unchanged by scaling q and b together, so store a unit q.
      const float norm = std::sqrt(dot4(quatAt(q), quatAt(q)));
      if (norm == 0.0f) {
        SpzLog("[SPZ ERROR] buildTemporalCloud: Gaussian %zu has a zero rotation at its origin", i);
        return std::nullopt;
      }
      for (size_t c = 0; c < kQuaternionComponents; ++c) q[c] /= norm;
      for (size_t k = 1; k <= source.rotationDegree; ++k)
        for (size_t c = 0; c < kQuaternionComponents; ++c)
          b[(k - 1) * kQuaternionComponents + c] *= static_cast<float>(spanPowers[k]) / norm;
    }
    bool moves = std::any_of(&posCoeffs[i * posDim], &posCoeffs[i * posDim] + posDim,
                             [](float v) { return v != 0.0f; }) ||
                 std::any_of(&rotCoeffs[i * rotDim], &rotCoeffs[i * rotDim] + rotDim,
                             [](float v) { return v != 0.0f; });
    if (row) {
      poster.alphas[i] = alphaLogit(windowOpacity(row, windowDim, source.posterTime));
      // Half the always-visible sigma: quantizing its log can't bring it below that.
      const bool fades = std::exp(static_cast<double>(row[kWindowLogSigma])) < alwaysVisibleSigma / 2;
      const bool clipped = hasInterval &&
                           (row[kWindowLo] > source.clipStart || row[kWindowHi] < source.clipEnd);
      moves = moves || fades || clipped;
    }
    animated[i] = moves;
  }

  // Animated Gaussians first, so channels cover a single leading range, in the order they become
  // visible, so each time bin's points are close together (see SpzTemporalTimeIndex).
  std::vector<double> visibleFrom(n, source.clipStart);
  if (windowed) {
    for (size_t i = 0; i < n; ++i) {
      const auto span = visibleSpanFor(&window[i * windowDim], windowDim, *clip);
      visibleFrom[i] = span ? span->first : std::numeric_limits<double>::infinity();
    }
  }
  std::vector<size_t> order;
  order.reserve(n);
  for (size_t i = 0; i < n; ++i)
    if (animated[i]) order.push_back(i);
  std::stable_sort(order.begin(), order.end(),
                   [&](size_t a, size_t b) { return visibleFrom[a] < visibleFrom[b]; });
  const size_t animatedCount = order.size();
  for (size_t i = 0; i < n; ++i)
    if (!animated[i]) order.push_back(i);

  const size_t shFloats = static_cast<size_t>(dimForDegree(g.shDegree)) * 3;
  GaussianCloud out = poster;
  for (size_t dst = 0; dst < n; ++dst) copyPoint(poster, order[dst], dst, shFloats, &out);

  const auto addChannel = [&](SpzTemporalTarget target, SpzTemporalChannelType type, size_t dim,
                              const std::vector<float>& values, std::vector<float> steps) {
    SpzTemporalChannel c;
    c.target = target;
    c.type = type;
    c.dim = static_cast<uint16_t>(dim);
    c.encoding = SpzChannelEncoding::Step;
    c.steps = std::move(steps);
    if (animatedCount < n) c.pointRange = SpzPointRange{0, static_cast<uint32_t>(animatedCount)};
    c.values.resize(animatedCount * dim);
    for (size_t dst = 0; dst < animatedCount; ++dst) {
      std::copy_n(&values[order[dst] * dim], dim, &c.values[dst * dim]);
    }
    clip->channels.push_back(std::move(c));
  };
  if (animatedCount > 0) {
    using T = SpzTemporalTarget;
    using Y = SpzTemporalChannelType;
    if (posDim > 0) {
      addChannel(T::Position, Y::Polynomial, posDim, posCoeffs,
                 std::vector<float>(posDim, positionStep));
    }
    if (rotDim > 0) {
      addChannel(T::Rotation, Y::Polynomial, rotDim, rotCoeffs,
                 std::vector<float>(rotDim, rotationStep));
    }
    if (windowed) addChannel(T::Opacity, Y::Window, windowDim, window, windowSteps);
  }
  out.extensions.push_back(clip);
  return out;
}

std::vector<float> channelMaxError(const SpzTemporalChannel& channel) {
  if (channel.dim == 0 || !encodingProblem(channel).empty()) return {};
  switch (channel.encoding) {
    case SpzChannelEncoding::Float32:
      return std::vector<float>(channel.dim, 0.0f);
    case SpzChannelEncoding::Range: {
      std::vector<float> min;
      std::vector<float> max;
      componentRanges(channel.values, channel.dim, &min, &max);
      const double levels = static_cast<double>((1u << channel.bits) - 1);
      std::vector<float> error(channel.dim);
      for (size_t c = 0; c < channel.dim; ++c) {
        error[c] = static_cast<float>((static_cast<double>(max[c]) - min[c]) / (2.0 * levels));
      }
      return error;
    }
    case SpzChannelEncoding::Step: {
      std::vector<float> error(channel.steps);
      for (float& e : error) e /= 2.0f;
      return error;
    }
  }
  return {};
}

}  // namespace spz
