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


#ifndef SPZ_TEMPORAL_NIANTIC_H_
#define SPZ_TEMPORAL_NIANTIC_H_

#include "splat-extensions.h"
#include "splat-types.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace spz {

// The property a channel animates.
enum class SpzTemporalTarget : uint8_t {
  Position = 1,
  Rotation = 2,
  Opacity = 3,
  // 4 (scale) and 5 (color) are reserved.
  Transform = 6,  // position and rotation together, driven by rigid drivers
};

// How a channel animates its target.
enum class SpzTemporalChannelType : uint8_t {
  Polynomial = 1,  // per-point polynomial in (t - posterTime)
  Window = 2,      // per-point Gaussian visibility window
  // 3 (per-point keyframes) is reserved.
  Skin = 4,  // per-point weighted blend of the record's rigid drivers
};

enum class SpzTemporalChannelFlags : uint8_t {
  None = 0,
  // Readers that don't support the channel must treat the cloud as static.
  Required = 1,
};

// Playback hint for a standalone clip; evaluation itself only takes clip time.
enum class SpzTemporalLoop : uint8_t { Once = 0, Loop = 1, PingPong = 2 };

// How a channel's values are stored.
enum class SpzChannelEncoding : uint8_t {
  Float32 = 0,  // lossless
  Range = 1,    // 1-16 bit levels between each component's min and max over the covered points
  Step = 2,     // integer multiples of a per-component step: error at most step / 2
};

// Default bits for Range channels.
constexpr uint8_t kDefaultRangeBits = 16;

// One animated property and its per-point data. Supported target/type pairs and their layouts
// (dim floats per covered point):
//
//   Position / Polynomial: dim 3·d (d = 1..3): [d1 xyz, d2 xyz, d3 xyz]
//   Rotation / Polynomial: dim 4·e (e = 1..3): xyzw terms added to the core rotation
//   Opacity / Window:      dim 3: [center (s), ln sigma (s), peak opacity in [0, 1]], or dim 5
//                          adding a hard visibility interval [lo, hi) (s)
//   Transform / Skin:      dim 2·J (J = 1..8): J driver indices, then J non-negative weights
//
// Polynomial terms multiply u^k, where u = (t - origin) / span is time normalized over the point's
// visible span (see temporalTimeFrame), so a term is the displacement it causes across that span.
// Values are point-major, for every point or only pointRange.
struct SpzTemporalChannel {
  SpzTemporalTarget target = SpzTemporalTarget::Position;
  SpzTemporalChannelType type = SpzTemporalChannelType::Polynomial;
  SpzTemporalChannelFlags flags = SpzTemporalChannelFlags::Required;
  uint16_t dim = 0;
  SpzChannelEncoding encoding = SpzChannelEncoding::Float32;
  // Range only: 1-16 (values must be finite).
  uint8_t bits = kDefaultRangeBits;
  // Step only: dim finite, positive steps (values must be finite and within 2^31 steps of 0).
  std::vector<float> steps;
  std::optional<SpzPointRange> pointRange;
  std::vector<float> values;
};

// The largest error saving can introduce into each of channel's components (0 for Float32), or an
// empty vector if channel's encoding parameters are invalid.
std::vector<float> channelMaxError(const SpzTemporalChannel& channel);

// Floats per driver keyframe: translation xyz, then rotation quaternion xyzw.
constexpr size_t kDriverTransformFloats = 7;

// Default opacityCutoff: one 8-bit alpha level.
constexpr float kDefaultOpacityCutoff = 1.0f / 255.0f;

// SPZ_NIANTIC_temporal: makes the cloud an animated clip (4D Gaussian splatting). The core
// attributes hold the scene at posterTime, so readers without this extension render a correct still
// frame; channels animate it over [clipStart, clipEnd). The record is self-contained: its channels
// and drivers keep their data in the record's own extension frames. See
// extensions/temporal-design.md for the evaluation equations.
//
// Skinning blends rigid drivers, each a keyframed transform from the bind pose (the core state)
// sampled at driverKeyframeTimes. A cloud has at most one temporal record, saved only in SPZ
// version 4 files.
//
// Preconditions for saving (checked; saveSpz fails and logs the problem otherwise):
//   - clipStart < clipEnd, clipStart <= posterTime < clipEnd, all finite.
//   - opacityCutoff in [0, 1); nominalFps finite and >= 0; loop a known mode.
//   - At most one channel per target; required channels use a supported target and type, with the
//     dim their layout needs.
//   - Each channel: a known encoding with valid parameters (Range: bits 1-16; Step: dim finite
//     positive steps); pointRange (if set) satisfies begin < end <= numPoints; values holds dim
//     floats per covered point (finite when quantized, within 2^31 steps for Step); skin indices in
//     [0, driverCount) and weights >= 0; window intervals finite.
//   - driverKeyframeTimes strictly increasing and finite; driverTransforms holds
//     driverCount * driverKeyframeTimes.size() * 7 finite floats with non-zero quaternions;
//     a Transform channel needs driverCount >= 1 and at least one keyframe; driverNames is empty or
//     holds driverCount names of 1-65535 printable ASCII bytes.
struct SpzExtensionTemporalNiantic : public SpzExtensionBase,
                                     public SpzFramedExtension,
                                     public SpzKeyedExtension,
                                     public SpzCoordinateDependentExtension,
                                     public SpzRotationCoupledExtension {
  double clipStart = 0.0;
  double clipEnd = 1.0;
  double posterTime = 0.0;
  SpzTemporalLoop loop = SpzTemporalLoop::Loop;
  float nominalFps = 0.0f;  // 0 = unspecified
  float opacityCutoff = kDefaultOpacityCutoff;
  std::vector<SpzTemporalChannel> channels;

  uint32_t driverCount = 0;
  std::vector<double> driverKeyframeTimes;
  // [driver][keyframe][tx ty tz qx qy qz qw]: the transform from the bind pose.
  std::vector<float> driverTransforms;
  // Optional labels (bone or object names), one per driver, for editing and retargeting.
  std::vector<std::string> driverNames;

  SpzExtensionTemporalNiantic();

  // The record's data lives in extension frames, which saveSpz writes through planFrameWrites().
  // write(os) therefore writes nothing (and logs an error).
  uint32_t payloadBytes() const override;
  void write(std::ostream& os) const override;
  SpzExtensionBase* copyAsRawData() const override;

  // PLY can't represent an animated clip: saving to PLY keeps the poster frame and drops the
  // record with a warning. Use evaluateTemporal to export a frame.
  std::optional<std::shared_ptr<SpzExtensionBase>> tryReadFromPly(
      std::istream& in, const std::unordered_set<std::string>& elementNames) const override;
  void writePlyHeader(std::ostream& out) const override;
  void writePlyData(std::ostream& out) const override;

  std::optional<SpzFrameWritePlan> planFrameWrites(int32_t numPoints) const override;
  std::string uniquenessKey() const override;
  // Converts drivers and position/rotation polynomial channels.
  SpzExtensionBasePtr withConvertedCoordinates(const CoordinateConverter& c) const override;
  // Rescales rotation polynomial terms with their point's stored core rotation.
  SpzExtensionBasePtr withRotationFactors(const std::vector<float>& factors) const override;

  // The number of keyframes per driver.
  size_t driverKeyframeCount() const { return driverKeyframeTimes.size(); }

  // The channel animating target, or nullptr.
  const SpzTemporalChannel* findChannel(SpzTemporalTarget target) const;

  // This extension's type ID, SPZ_NIANTIC_temporal.
  static SpzExtensionType type();
};

// A temporal record validated against its cloud, with its supported channels looked up. Absent
// channels are nullptr; the pointers live as long as clip.
struct SpzTemporalModel {
  std::shared_ptr<const SpzExtensionTemporalNiantic> clip;
  const SpzTemporalChannel* position = nullptr;
  const SpzTemporalChannel* rotation = nullptr;
  const SpzTemporalChannel* opacity = nullptr;
  const SpzTemporalChannel* skin = nullptr;
};

// Binds the cloud's temporal record. Returns std::nullopt when there is none, or (after logging
// why) when it doesn't satisfy the saving preconditions for this cloud, which must then be treated
// as static.
std::optional<SpzTemporalModel> bindTemporalModel(const GaussianCloud& cloud);

// Where a point's polynomials are measured from, and over what span. With the point's visible span
// [lo, hi] — the clip, intersected with its window's support |t - center| <= K·sigma (K =
// sqrt(-2 ln opacityCutoff)) and its hard interval when it has them — origin is posterTime clamped
// into [lo, hi] and span is hi - lo. A point visible at posterTime therefore has origin posterTime,
// and the core attributes stay the poster frame. A point never visible gets posterTime and the
// clip's duration.
struct SpzTemporalTimeFrame {
  double origin = 0.0;
  double span = 1.0;
};

// The time frame of point under model. Precondition: model is bindTemporalModel of a cloud with
// more than point points.
SpzTemporalTimeFrame temporalTimeFrame(const SpzTemporalModel& model, size_t point);

// Which points can be visible when: the clip divided into time bins [boundaries[b],
// boundaries[b + 1]), each with the range of points that may be at least opacityCutoff opaque
// during it (a superset: points outside a bin's range are invisible throughout it). Writers sort
// points by when they become visible, so ranges are tight; a renderer then evaluates only the bin
// holding t.
struct SpzTemporalTimeIndex {
  std::vector<double> boundaries;  // bins.size() + 1 increasing times, from clipStart to clipEnd
  std::vector<SpzPointRange> bins;
  uint32_t maxBinPoints = 0;       // the largest bin's point count: enough slots for any frame
};

// The bin boundaries writers use for model's clip: its windows' distinct interval bounds when it
// has intervals (so a frame sequence gets exactly one bin per frame), else equal bins, one per frame
// at its nominalFps or a default count.
std::vector<double> temporalTimeBinBoundaries(const SpzTemporalModel& model);

// count + 1 equally spaced boundaries over [start, end). Precondition: count >= 1, start < end.
std::vector<double> uniformTimeBinBoundaries(double start, double end, size_t count);

// Computes the time index of the cloud model was bound to over boundaries. Preconditions: model is
// bindTemporalModel of a cloud of numPoints points; boundaries holds at least 2 increasing times
// from clipStart to clipEnd.
SpzTemporalTimeIndex computeTemporalTimeIndex(const SpzTemporalModel& model, int32_t numPoints,
                                              const std::vector<double>& boundaries);

// The points that may be visible at clip time t (empty outside the clip).
SpzPointRange temporalCandidates(const SpzTemporalTimeIndex& index, double t);

// The attributes that change over time, for one frame. Reuse one across frames so evaluating
// doesn't allocate.
struct SpzTemporalFrame {
  std::vector<float> positions;  // numPoints * 3
  std::vector<float> rotations;  // numPoints * 4, xyzw
  std::vector<float> alphas;     // numPoints, as logits (like GaussianCloud::alphas)
};

// Writes the positions, rotations and alphas at clip time t into *frame; scales, colors and SH
// don't change over time, so callers keep using the cloud's. Precondition: model is
// bindTemporalModel(cloud).
void evaluateTemporalInto(const GaussianCloud& cloud, const SpzTemporalModel& model, double t,
                          SpzTemporalFrame* frame);

// The Gaussians visible (at least opacityCutoff opaque) at one time, packed in point order. Reuse
// one across frames so evaluating doesn't allocate.
struct SpzTemporalVisibleFrame {
  std::vector<uint32_t> indices;  // the visible points, ascending
  std::vector<float> positions;   // indices.size() * 3
  std::vector<float> rotations;   // indices.size() * 4, xyzw
  std::vector<float> alphas;      // indices.size(), as logits
  SpzTemporalFrame work;          // working storage, reused across frames
};

// Evaluates only the points index says may be visible at clip time t, and keeps those that are.
// Scales, colors and SH don't change: gather them through frame->indices. Preconditions: model is
// bindTemporalModel(cloud); index is computeTemporalTimeIndex of model and cloud.numPoints.
void evaluateVisibleTemporal(const GaussianCloud& cloud, const SpzTemporalModel& model,
                             const SpzTemporalTimeIndex& index, double t,
                             SpzTemporalVisibleFrame* frame);

// Returns the static cloud at clip time t, without the temporal record (other extensions are
// kept). Gaussians outside [clipStart, clipEnd) are fully transparent. Precondition: model is
// bindTemporalModel(cloud).
GaussianCloud evaluateTemporal(const GaussianCloud& cloud, const SpzTemporalModel& model, double t);

// Maps seconds of playback to clip time using the clip's loop mode. Precondition: clip satisfies
// the record's save preconditions (clipStart < clipEnd).
double clipTime(const SpzExtensionTemporalNiantic& clip, double playbackSeconds);

// The largest error the builder may introduce into each stored quantity (it uses Step channels with
// steps of twice these bounds). A polynomial's total error is at most the sum of its terms' bounds.
struct SpzTemporalPrecision {
  float position = 5e-4f;        // per position polynomial term, scene units
  float rotation = 1e-3f;        // per rotation polynomial term component
  float time = 5e-5f;            // window centers and interval bounds, seconds
  float logSigma = 1e-3f;        // ln of window sigma (about 0.1% of sigma)
  float opacity = 1.0f / 512.0f; // window peak
};

// Per-Gaussian motion in the convention most 4DGS methods train in: each Gaussian is described at
// its own reference time tCenter, with polynomials in (t - tCenter) and a Gaussian visibility
// window centered at tCenter.
struct SpzTemporalSource {
  // Each Gaussian's state at its tCenter; alphas hold the window's peak opacity (as a logit).
  // Must not carry custom attributes (points are reordered) or a temporal record.
  GaussianCloud cloud;
  // numPoints reference times in seconds; empty means every Gaussian's is posterTime.
  std::vector<float> tCenter;
  // numPoints ln(sigma) of the visibility window in seconds; empty means always visible.
  std::vector<float> tLogSigma;
  // numPoints * 2 hard visibility intervals [lo, hi) in seconds; empty means none.
  std::vector<float> tInterval;
  // 0-3; numPoints * 3 * positionDegree floats: [a1 xyz, a2 xyz, a3 xyz] per point.
  uint8_t positionDegree = 0;
  std::vector<float> positionCoefficients;
  // 0-3; numPoints * 4 * rotationDegree floats: xyzw terms per degree per point.
  uint8_t rotationDegree = 0;
  std::vector<float> rotationCoefficients;

  double clipStart = 0.0;
  double clipEnd = 1.0;
  double posterTime = 0.0;
  SpzTemporalLoop loop = SpzTemporalLoop::Loop;
  float nominalFps = 0.0f;
  float opacityCutoff = kDefaultOpacityCutoff;
  SpzTemporalPrecision precision;
};

// Builds an animated cloud from source: re-expresses the polynomials in each point's time frame,
// writes each point's state at its origin into the core attributes (the poster frame for every
// point visible then), moves animated Gaussians before static ones (channels then cover only the
// animated range), and adds the temporal record with Step channels at source.precision. Returns std::nullopt (after
// logging why) if source is inconsistent.
std::optional<GaussianCloud> buildTemporalCloud(const SpzTemporalSource& source);

}  // namespace spz

#endif  // SPZ_TEMPORAL_NIANTIC_H_
