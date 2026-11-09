// Extension bindings for the spz Python module.
// This file contains all Python bindings related to SPZ extensions.

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>

#include <cstring>

#include "src/cc/load-spz.h"
#include "src/cc/splat-types.h"
#include "extensions/cc/splat-extensions.h"
#include "extensions/cc/safe-orbit-camera-adobe.h"
#include "extensions/cc/coordinate-system-adobe.h"
#include "extensions/cc/georeference-niantic.h"
#include "extensions/cc/custom-attribute-niantic.h"
#include "extensions/cc/temporal-niantic.h"
#include "extensions/python/splat-extensions.h"

namespace nb = nanobind;

namespace {
// Any-shape, C-contiguous float32 array, so values can be set from e.g. a (num_points, dim) array.
using FloatArray = nb::ndarray<nb::numpy, float, nb::c_contig, nb::device::cpu>;

// Returns a NumPy array that owns a copy of values.
template <typename T>
nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig, nb::device::cpu> copyToNumpy(
    const std::vector<T>& values) {
    using Array = nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    if (values.empty()) return Array(nullptr, {0}, nb::handle());
    T* data = new T[values.size()];
    std::memcpy(data, values.data(), values.size() * sizeof(T));
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<T*>(p); });
    return Array(data, {values.size()}, owner);
}
}  // namespace

namespace spz {
namespace python {
namespace {
// Copies an optional NumPy array into a vector, or leaves it empty for None.
std::vector<float> toVector(const std::optional<FloatArray>& arr) {
    if (!arr) return {};
    return std::vector<float>(arr->data(), arr->data() + arr->size());
}

void register_temporal(nb::module_& m) {
    nb::enum_<spz::SpzTemporalTarget>(m, "SpzTemporalTarget", "The property a temporal channel animates.")
        .value("POSITION", spz::SpzTemporalTarget::Position)
        .value("ROTATION", spz::SpzTemporalTarget::Rotation)
        .value("OPACITY", spz::SpzTemporalTarget::Opacity)
        .value("TRANSFORM", spz::SpzTemporalTarget::Transform,
               "Position and rotation together, driven by the record's rigid drivers");
    nb::enum_<spz::SpzTemporalChannelType>(m, "SpzTemporalChannelType", "How a temporal channel animates its target.")
        .value("POLYNOMIAL", spz::SpzTemporalChannelType::Polynomial, "Per-point polynomial in (t - poster_time)")
        .value("WINDOW", spz::SpzTemporalChannelType::Window, "Per-point Gaussian visibility window")
        .value("SKIN", spz::SpzTemporalChannelType::Skin, "Per-point weighted blend of rigid drivers");
    nb::enum_<spz::SpzTemporalChannelFlags>(m, "SpzTemporalChannelFlags")
        .value("NONE", spz::SpzTemporalChannelFlags::None)
        .value("REQUIRED", spz::SpzTemporalChannelFlags::Required,
               "Readers that don't support the channel must treat the cloud as static");
    nb::enum_<spz::SpzTemporalLoop>(m, "SpzTemporalLoop", "Playback hint for a standalone clip.")
        .value("ONCE", spz::SpzTemporalLoop::Once)
        .value("LOOP", spz::SpzTemporalLoop::Loop)
        .value("PING_PONG", spz::SpzTemporalLoop::PingPong);

    nb::enum_<spz::SpzChannelEncoding>(m, "SpzChannelEncoding", "How a temporal channel's values are stored.")
        .value("FLOAT32", spz::SpzChannelEncoding::Float32, "Lossless")
        .value("RANGE", spz::SpzChannelEncoding::Range,
               "bits-bit levels between each component's min and max over the covered points")
        .value("STEP", spz::SpzChannelEncoding::Step,
               "Integer multiples of a per-component step: error at most step / 2");

    nb::class_<spz::SpzTemporalChannel>(m, "SpzTemporalChannel",
        R"doc(One animated property and its per-point data (dim floats per covered point, point-major).

        POSITION/POLYNOMIAL: dim 3·d, [d1 xyz, d2 xyz, d3 xyz]. ROTATION/POLYNOMIAL: dim 4·e, xyzw terms
        added to the core rotation. Polynomial terms multiply u^k, with u the time normalized over the
        point's visible span (see temporal_time_frame). OPACITY/WINDOW: dim 3, [center s, ln sigma s,
        peak opacity], or dim 5 adding a hard interval [lo, hi) s. TRANSFORM/SKIN: dim 2·J, J driver
        indices then J weights. encoding is FLOAT32 (default), RANGE (bits 1-16) or STEP (steps, one per
        component); point_range limits the channel to points begin..end-1.)doc")
        .def(nb::init<>())
        .def("__init__",
             [](spz::SpzTemporalChannel* self, spz::SpzTemporalTarget target,
                spz::SpzTemporalChannelType type, uint16_t dim, const FloatArray& values,
                spz::SpzChannelEncoding encoding, uint8_t bits, const std::optional<FloatArray>& steps,
                const std::optional<std::tuple<uint32_t, uint32_t>>& point_range,
                spz::SpzTemporalChannelFlags flags) {
                 new (self) spz::SpzTemporalChannel();
                 self->target = target;
                 self->type = type;
                 self->dim = dim;
                 self->values.assign(values.data(), values.data() + values.size());
                 self->encoding = encoding;
                 self->bits = bits;
                 self->steps = toVector(steps);
                 if (point_range) {
                     self->pointRange = spz::SpzPointRange{std::get<0>(*point_range), std::get<1>(*point_range)};
                 }
                 self->flags = flags;
             },
             nb::arg("target"), nb::arg("type"), nb::arg("dim"), nb::arg("values"),
             nb::arg("encoding") = spz::SpzChannelEncoding::Float32,
             nb::arg("bits") = spz::kDefaultRangeBits, nb::arg("steps").none() = nb::none(),
             nb::arg("point_range").none() = nb::none(),
             nb::arg("flags") = spz::SpzTemporalChannelFlags::Required)
        .def_rw("target", &spz::SpzTemporalChannel::target)
        .def_rw("type", &spz::SpzTemporalChannel::type)
        .def_rw("flags", &spz::SpzTemporalChannel::flags)
        .def_rw("dim", &spz::SpzTemporalChannel::dim)
        .def_rw("encoding", &spz::SpzTemporalChannel::encoding)
        .def_rw("bits", &spz::SpzTemporalChannel::bits, "RANGE only: 1-16")
        .def_prop_rw("steps",
                     [](const spz::SpzTemporalChannel& self) { return copyToNumpy(self.steps); },
                     [](spz::SpzTemporalChannel& self, const FloatArray& arr) {
                         self.steps.assign(arr.data(), arr.data() + arr.size());
                     },
                     nb::rv_policy::move, "STEP only: one positive step per component (copied on access)")
        .def_prop_rw("point_range",
                     [](const spz::SpzTemporalChannel& self) -> std::optional<std::tuple<uint32_t, uint32_t>> {
                         if (!self.pointRange) return std::nullopt;
                         return std::make_tuple(self.pointRange->begin, self.pointRange->end);
                     },
                     [](spz::SpzTemporalChannel& self, const std::optional<std::tuple<uint32_t, uint32_t>>& range) {
                         if (!range) {
                             self.pointRange.reset();
                             return;
                         }
                         self.pointRange = spz::SpzPointRange{std::get<0>(*range), std::get<1>(*range)};
                     },
                     nb::arg("range").none(), "(begin, end), or None for every point")
        .def_prop_rw("values",
                     [](const spz::SpzTemporalChannel& self) { return copyToNumpy(self.values); },
                     [](spz::SpzTemporalChannel& self, const FloatArray& arr) {
                         self.values.assign(arr.data(), arr.data() + arr.size());
                     },
                     nb::rv_policy::move, "dim floats per covered point, point-major (copied on access)");

    m.def("channel_max_error",
          [](const spz::SpzTemporalChannel& channel) { return copyToNumpy(spz::channelMaxError(channel)); },
          nb::arg("channel"),
          "The largest error saving can introduce into each of the channel's components (0 for "
          "FLOAT32); empty if its encoding parameters are invalid.");

    nb::class_<spz::SpzExtensionTemporalNiantic, spz::SpzExtensionBase>(m, "SpzExtensionTemporalNiantic",
        R"doc(Makes the cloud an animated clip (4D Gaussian splatting), saved only in SPZ v4 files.

        The core attributes hold the scene at poster_time; channels animate it over
        [clip_start, clip_end), each holding its own per-point data.
        Skinned clips (a TRANSFORM/SKIN channel) blend driver_count rigid drivers, keyframed at
        driver_keyframe_times; driver_transforms holds driver_count * keyframes * 7 floats
        (translation xyz, rotation quaternion xyzw) transforming the bind pose (the core state).
        See extensions/temporal-design.md for the evaluation equations.)doc")
        .def(nb::init<>())
        .def_rw("clip_start", &spz::SpzExtensionTemporalNiantic::clipStart, "Clip start, seconds")
        .def_rw("clip_end", &spz::SpzExtensionTemporalNiantic::clipEnd, "Clip end (exclusive), seconds")
        .def_rw("poster_time", &spz::SpzExtensionTemporalNiantic::posterTime,
                "Time the core attributes show, in [clip_start, clip_end)")
        .def_rw("loop", &spz::SpzExtensionTemporalNiantic::loop, "Playback hint")
        .def_rw("nominal_fps", &spz::SpzExtensionTemporalNiantic::nominalFps,
                "Playback hint, 0 = unspecified; evaluation is continuous")
        .def_rw("opacity_cutoff", &spz::SpzExtensionTemporalNiantic::opacityCutoff,
                "Renderers may skip Gaussians whose opacity at t is below this")
        .def_rw("channels", &spz::SpzExtensionTemporalNiantic::channels,
                "Animated properties; assign a new list, in-place edits don't stick")
        .def_rw("driver_count", &spz::SpzExtensionTemporalNiantic::driverCount, "Number of rigid drivers")
        .def_rw("driver_keyframe_times", &spz::SpzExtensionTemporalNiantic::driverKeyframeTimes,
                "Strictly increasing keyframe times shared by every driver, seconds")
        .def_rw("driver_names", &spz::SpzExtensionTemporalNiantic::driverNames,
                "Optional labels, one per driver (bone or object names); assign a new list")
        .def_prop_rw("driver_transforms",
                     [](const spz::SpzExtensionTemporalNiantic& self) { return copyToNumpy(self.driverTransforms); },
                     [](spz::SpzExtensionTemporalNiantic& self, const FloatArray& arr) {
                         self.driverTransforms.assign(arr.data(), arr.data() + arr.size());
                     },
                     nb::rv_policy::move,
                     "driver_count * keyframes * 7 floats [tx ty tz qx qy qz qw], driver-major; set from "
                     "any C-contiguous array, e.g. (driver_count, keyframes, 7) (copied on access)")
        .def_static("type", &spz::SpzExtensionTemporalNiantic::type);

    m.def("find_temporal",
          [](const spz::GaussianCloud& cloud) {
              return spz::findExtensionByType<spz::SpzExtensionTemporalNiantic>(cloud.extensions);
          },
          nb::arg("cloud"), "Returns the cloud's temporal record, or None.");
    m.def("is_animated",
          [](const spz::GaussianCloud& cloud) { return spz::bindTemporalModel(cloud).has_value(); },
          nb::arg("cloud"),
          "True if the cloud has a valid temporal record (otherwise the cloud is static; the "
          "reason is logged).");
    m.def("evaluate_temporal",
          [](const spz::GaussianCloud& cloud, double t) {
              std::optional<spz::SpzTemporalModel> model = spz::bindTemporalModel(cloud);
              if (!model) throw std::invalid_argument("the cloud isn't a valid animated clip");
              return spz::evaluateTemporal(cloud, *model, t);
          },
          nb::arg("cloud"), nb::arg("t"),
          "Returns the static cloud at clip time t (seconds), without the temporal record. Raises "
          "ValueError if the cloud isn't a valid animated clip.");
    m.def("clip_time", &spz::clipTime, nb::arg("clip"), nb::arg("playback_seconds"),
          "Maps seconds of playback to clip time using the clip's loop mode.");
    m.def("temporal_time_index",
          [](const spz::GaussianCloud& cloud, std::optional<size_t> bins) {
              std::optional<spz::SpzTemporalModel> model = spz::bindTemporalModel(cloud);
              if (!model) throw std::invalid_argument("the cloud isn't a valid animated clip");
              if (bins && *bins == 0) throw std::invalid_argument("bins must be at least 1");
              const std::vector<double> boundaries =
                  bins ? spz::uniformTimeBinBoundaries(model->clip->clipStart, model->clip->clipEnd, *bins)
                       : spz::temporalTimeBinBoundaries(*model);
              const spz::SpzTemporalTimeIndex index =
                  spz::computeTemporalTimeIndex(*model, cloud.numPoints, boundaries);
              std::vector<std::tuple<uint32_t, uint32_t>> ranges;
              for (const spz::SpzPointRange& r : index.bins) ranges.emplace_back(r.begin, r.end);
              return std::make_tuple(index.boundaries, ranges, index.maxBinPoints);
          },
          nb::arg("cloud"), nb::arg("bins").none() = nb::none(),
          "(boundaries, [(begin, end) per bin], max_bin_points): the points that may be visible in "
          "each time bin [boundaries[b], boundaries[b + 1]). bins=None uses the boundaries writers "
          "store (interval bounds, else one bin per frame at the nominal fps, else 64); an integer "
          "uses that many equal bins.");
    m.def("evaluate_visible_temporal",
          [](const spz::GaussianCloud& cloud, double t) {
              std::optional<spz::SpzTemporalModel> model = spz::bindTemporalModel(cloud);
              if (!model) throw std::invalid_argument("the cloud isn't a valid animated clip");
              const spz::SpzTemporalTimeIndex index = spz::computeTemporalTimeIndex(
                  *model, cloud.numPoints, spz::temporalTimeBinBoundaries(*model));
              spz::SpzTemporalVisibleFrame frame;
              spz::evaluateVisibleTemporal(cloud, *model, index, t, &frame);
              return std::make_tuple(copyToNumpy(frame.indices), copyToNumpy(frame.positions),
                                     copyToNumpy(frame.rotations), copyToNumpy(frame.alphas));
          },
          nb::arg("cloud"), nb::arg("t"),
          "(indices, positions, rotations, alphas) of the Gaussians at least opacity_cutoff opaque at "
          "clip time t, evaluating only the points the time index says may be.");
    m.def("temporal_time_frame",
          [](const spz::GaussianCloud& cloud, size_t point) {
              std::optional<spz::SpzTemporalModel> model = spz::bindTemporalModel(cloud);
              if (!model) throw std::invalid_argument("the cloud isn't a valid animated clip");
              if (point >= static_cast<size_t>(cloud.numPoints)) throw std::out_of_range("point");
              const spz::SpzTemporalTimeFrame frame = spz::temporalTimeFrame(*model, point);
              return std::make_tuple(frame.origin, frame.span);
          },
          nb::arg("cloud"), nb::arg("point"),
          "(origin, span) of a point's polynomials: u = (t - origin) / span. origin is poster_time "
          "for every point visible then.");
    m.def("build_temporal_cloud",
          [](const spz::GaussianCloud& cloud, double clip_start, double clip_end, double poster_time,
             const std::optional<FloatArray>& t_center, const std::optional<FloatArray>& t_log_sigma,
             const std::optional<FloatArray>& t_interval,
             const std::optional<FloatArray>& position_coefficients, uint8_t position_degree,
             const std::optional<FloatArray>& rotation_coefficients, uint8_t rotation_degree,
             spz::SpzTemporalLoop loop, float nominal_fps, float opacity_cutoff,
             float position_error, float rotation_error, float time_error, float log_sigma_error,
             float opacity_error) {
              spz::SpzTemporalSource source;
              source.cloud = cloud;
              source.tCenter = toVector(t_center);
              source.tLogSigma = toVector(t_log_sigma);
              source.tInterval = toVector(t_interval);
              source.positionCoefficients = toVector(position_coefficients);
              source.positionDegree = position_degree;
              source.rotationCoefficients = toVector(rotation_coefficients);
              source.rotationDegree = rotation_degree;
              source.clipStart = clip_start;
              source.clipEnd = clip_end;
              source.posterTime = poster_time;
              source.loop = loop;
              source.nominalFps = nominal_fps;
              source.opacityCutoff = opacity_cutoff;
              source.precision = {position_error, rotation_error, time_error, log_sigma_error,
                                  opacity_error};
              std::optional<spz::GaussianCloud> built = spz::buildTemporalCloud(source);
              if (!built) throw std::invalid_argument("inconsistent temporal source (see log)");
              return std::move(*built);
          },
          nb::arg("cloud"), nb::arg("clip_start"), nb::arg("clip_end"), nb::arg("poster_time"),
          nb::arg("t_center").none() = nb::none(), nb::arg("t_log_sigma").none() = nb::none(),
          nb::arg("t_interval").none() = nb::none(),
          nb::arg("position_coefficients").none() = nb::none(), nb::arg("position_degree") = 0,
          nb::arg("rotation_coefficients").none() = nb::none(), nb::arg("rotation_degree") = 0,
          nb::arg("loop") = spz::SpzTemporalLoop::Loop, nb::arg("nominal_fps") = 0.0f,
          nb::arg("opacity_cutoff") = spz::kDefaultOpacityCutoff,
          nb::arg("position_error") = spz::SpzTemporalPrecision{}.position,
          nb::arg("rotation_error") = spz::SpzTemporalPrecision{}.rotation,
          nb::arg("time_error") = spz::SpzTemporalPrecision{}.time,
          nb::arg("log_sigma_error") = spz::SpzTemporalPrecision{}.logSigma,
          nb::arg("opacity_error") = spz::SpzTemporalPrecision{}.opacity,
          R"doc(Builds an animated clip from per-Gaussian motion given at each Gaussian's own time.

          cloud holds each Gaussian's state at its t_center (alphas: the window's peak opacity, as a
          logit) and must not carry custom attributes (points are reordered). t_center: (N,) seconds,
          or None for poster_time. t_log_sigma: (N,) ln(sigma) of the visibility window, or None for
          always visible. t_interval: (N, 2) hard visibility intervals [lo, hi), or None.
          position_coefficients: (N, 3 * position_degree) polynomial terms in (t - t_center),
          [a1 xyz, a2 xyz, a3 xyz]; rotation_coefficients likewise with 4 (xyzw) per degree. The
          *_error arguments bound the error saving may add to each stored quantity (position and
          rotation per polynomial term). Raises ValueError if inputs are inconsistent.)doc");
}
}  // namespace

// Register all extension-related Python bindings
void register_extensions(nb::module_& m) {
    // -------------------------------------------------------------------------
    // Extensions structs
    // -------------------------------------------------------------------------
    nb::enum_<spz::SpzExtensionType>(m, "SpzExtensionType", R"doc(
        Enumeration of vendor-specific SPZ extension types.
        More extensions may be added in the future, with their values defined as
        VENDOR_ID << 16 | EXTENSION_ID, where:
        - VENDOR_ID is a unique identifier for the vendor (e.g., Adobe = 0xADBE);
        - EXTENSION_ID is a unique identifier for the specific extension within that vendor's namespace.
    )doc")
        .value("SPZ_ADOBE_safe_orbit_camera", spz::SpzExtensionType::SPZ_ADOBE_safe_orbit_camera,
               "Adobe safe orbit camera extension")
        .value("SPZ_ADOBE_coordinate_system", spz::SpzExtensionType::SPZ_ADOBE_coordinate_system,
               "Adobe coordinate system extension — records the native coordinate system of the asset")
        .value("SPZ_NIANTIC_georeference", spz::SpzExtensionType::SPZ_NIANTIC_georeference,
               "Niantic georeference extension — similarity transform from the local frame to a body-fixed geocentric CRS")
        .value("SPZ_NIANTIC_custom_attribute", spz::SpzExtensionType::SPZ_NIANTIC_custom_attribute,
               "Niantic custom attribute extension — a named per-point attribute")
        .value("SPZ_NIANTIC_temporal", spz::SpzExtensionType::SPZ_NIANTIC_temporal,
               "Niantic temporal extension — an animated clip (4D Gaussian splatting)")
        .export_values();

    nb::class_<spz::SpzExtensionBase>(m, "SpzExtensionBase")
        .def_ro("extension_type", &spz::SpzExtensionBase::extensionType,
                "Type of the SPZ extension");

    nb::class_<spz::SpzExtensionSafeOrbitCameraAdobe, spz::SpzExtensionBase>(m, "SpzExtensionSafeOrbitCameraAdobe")
        .def(nb::init<>())
        .def_rw("safe_orbit_elevation_min", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMin,
                "Minimum elevation angle for safe orbit (radians)")
        .def_rw("safe_orbit_elevation_max", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMax,
                "Maximum elevation angle for safe orbit (radians)")
        .def_rw("safe_orbit_radius_min", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitRadiusMin,
                "Minimum radius for safe orbit")
        .def_static("type", &spz::SpzExtensionSafeOrbitCameraAdobe::type,
                    "Static method to get the extension type enum value");
    nb::class_<spz::SpzExtensionCoordinateSystemAdobe, spz::SpzExtensionBase>(m, "SpzExtensionCoordinateSystemAdobe")
        .def(nb::init<>())
        .def_rw("coordinate_system", &spz::SpzExtensionCoordinateSystemAdobe::coordinateSystem,
                "Native coordinate system of the asset's Gaussian data")
        .def_static("type", &spz::SpzExtensionCoordinateSystemAdobe::type,
                    "Static method to get the extension type enum value");
    nb::enum_<spz::SpzExtensionGeoreferenceNiantic::CrsEncoding>(
        m, "GeoreferenceCrsEncoding", "Serialization of SpzExtensionGeoreferenceNiantic.crs")
        .value("AUTHORITY_CODE", spz::SpzExtensionGeoreferenceNiantic::CrsEncoding::AuthorityCode)
        .value("PROJJSON", spz::SpzExtensionGeoreferenceNiantic::CrsEncoding::Projjson);
    nb::class_<spz::SpzExtensionGeoreferenceNiantic, spz::SpzExtensionBase>(m, "SpzExtensionGeoreferenceNiantic")
        .def(nb::init<>())
        .def_rw("crs", &spz::SpzExtensionGeoreferenceNiantic::crs,
                "Target geocentric CRS, e.g. \"EPSG:4978\" or inline PROJJSON; required")
        .def_rw("crs_encoding", &spz::SpzExtensionGeoreferenceNiantic::crsEncoding,
                "Serialization of crs: AUTHORITY_CODE (default) or PROJJSON (only when no code exists)")
        // Tuples, so element assignment raises instead of silently editing a copy.
        .def_prop_rw("origin",
                     [](const spz::SpzExtensionGeoreferenceNiantic& e) {
                       return std::make_tuple(e.origin[0], e.origin[1], e.origin[2]);
                     },
                     [](spz::SpzExtensionGeoreferenceNiantic& e, const std::array<double, 3>& v) { e.origin = v; },
                     "Translation from the RUB origin to the target CRS, in meters (x, y, z)")
        .def_prop_rw("rotation",
                     [](const spz::SpzExtensionGeoreferenceNiantic& e) {
                       return std::make_tuple(e.rotation[0], e.rotation[1], e.rotation[2], e.rotation[3]);
                     },
                     [](spz::SpzExtensionGeoreferenceNiantic& e, const std::array<double, 4>& v) { e.rotation = v; },
                     "Unit Hamilton quaternion (x, y, z, w) rotating RUB positions into the target CRS; "
                     "near-unit values are normalized on save")
        .def_rw("scale", &spz::SpzExtensionGeoreferenceNiantic::scale,
                "Dimensionless uniform scale from RUB positions into the target CRS; must be finite and > 0")
        .def_rw("epoch", &spz::SpzExtensionGeoreferenceNiantic::epoch,
                "Coordinate epoch of the transformed coordinates, as a decimal year; NaN when absent")
        .def_static("type", &spz::SpzExtensionGeoreferenceNiantic::type,
                    "Static method to get the extension type enum value");
    nb::enum_<spz::SpzAttributeSemantic>(m, "SpzAttributeSemantic",
        "How a custom attribute responds to coordinate conversions.")
        .value("CUSTOM", spz::SpzAttributeSemantic::Custom, "Never converted")
        .value("SPATIAL", spz::SpzAttributeSemantic::Spatial,
               "xyz triples in the cloud's coordinate space (normals, directions, ...), converted like positions");
    nb::class_<spz::SpzExtensionCustomAttributeNiantic, spz::SpzExtensionBase>(m, "SpzExtensionCustomAttributeNiantic",
        R"doc(A named per-point attribute, stored in SPZ v4 files in extension frames (zstd skippable frames).

        dim is the number of floats per point and semantic how they convert: SPATIAL needs dim to
        be a multiple of 3 and each xyz triple is converted with the cloud's coordinates; CUSTOM is
        never converted. point_range limits the attribute to a range of points. name is a free-form label. Saving fails unless names are 1-65535 bytes of printable
        ASCII and unique in the cloud, bits is 32 (lossless float32) or 1-16 (quantized per component), values
        holds num_points * dim floats, point-major, and the file version is 4 or later.)doc")
        .def(nb::init<>())
        .def_rw("name", &spz::SpzExtensionCustomAttributeNiantic::name, "Free-form label of printable ASCII, unique per cloud")
        .def_rw("semantic", &spz::SpzExtensionCustomAttributeNiantic::semantic,
                "SPATIAL values are converted with the cloud's coordinates; CUSTOM never are")
        .def_rw("dim", &spz::SpzExtensionCustomAttributeNiantic::dim, "Floats per point (1-65535)")
        .def_rw("bits", &spz::SpzExtensionCustomAttributeNiantic::bits,
                "32 for float32, or 1-16 for quantized storage")
        .def_prop_rw("point_range",
                     [](const spz::SpzExtensionCustomAttributeNiantic& self)
                         -> std::optional<std::tuple<uint32_t, uint32_t>> {
                         if (!self.pointRange) return std::nullopt;
                         return std::make_tuple(self.pointRange->begin, self.pointRange->end);
                     },
                     [](spz::SpzExtensionCustomAttributeNiantic& self,
                        const std::optional<std::tuple<uint32_t, uint32_t>>& range) {
                         if (!range) {
                             self.pointRange.reset();
                             return;
                         }
                         self.pointRange = spz::SpzPointRange{std::get<0>(*range), std::get<1>(*range)};
                     },
                     nb::arg("range").none(),
                     "(begin, end): the attribute covers only points begin..end-1, and values holds "
                     "(end - begin) * dim floats; None (default) covers every point")
        .def_prop_rw("values",
                     [](const spz::SpzExtensionCustomAttributeNiantic& self) { return copyToNumpy(self.values); },
                     [](spz::SpzExtensionCustomAttributeNiantic& self, const FloatArray& arr) {
                         self.values.assign(arr.data(), arr.data() + arr.size());
                     },
                     nb::rv_policy::move, "num_points * dim float32 values, point-major; set from any C-contiguous array, e.g. (num_points, dim) (copied on access)")
        .def_static("type", &spz::SpzExtensionCustomAttributeNiantic::type,
                    "Static method to get the extension type enum value");
    m.def("find_custom_attribute",
          [](const spz::GaussianCloud& cloud, const std::string& name) {
              return spz::findCustomAttribute(cloud.extensions, name);
          },
          nb::arg("cloud"), nb::arg("name"),
          "Returns the cloud's custom attribute called name, or None. Coordinate conversions "
          "replace SPATIAL attributes with converted copies, so look it up again after converting.");
    register_temporal(m);
    m.def("is_known_ply_extension_element", &spz::isKnownPlyExtensionElement, nb::arg("element_name"),
          "Returns True if the PLY extra element name is handled by an extension.");
}

// Register extension-related properties for GaussianCloud
void register_gaussian_cloud_extensions(nb::class_<spz::GaussianCloud>& cloud) {
    cloud.def_rw("extensions", &spz::GaussianCloud::extensions,
                 "List of SPZ extensions associated with this GaussianCloud.");
}

}  // namespace python
}  // namespace spz

