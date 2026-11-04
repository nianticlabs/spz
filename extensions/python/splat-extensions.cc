// Extension bindings for the spz Python module.
// This file contains all Python bindings related to SPZ extensions.

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
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
#include "extensions/cc/environment-map-adobe.h"
#include "extensions/python/splat-extensions.h"

namespace nb = nanobind;

namespace {
// Any-shape, C-contiguous float32 array, so values can be set from e.g. a (num_points, dim) array.
using FloatArray = nb::ndarray<nb::numpy, float, nb::c_contig, nb::device::cpu>;
using FloatArray1D = nb::ndarray<nb::numpy, float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// Returns a NumPy array that owns a copy of values.
FloatArray1D copyToNumpy(const std::vector<float>& values) {
    if (values.empty()) return FloatArray1D(nullptr, {0}, nb::handle());
    float* data = new float[values.size()];
    std::memcpy(data, values.data(), values.size() * sizeof(float));
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
    return FloatArray1D(data, {values.size()}, owner);
}
}  // namespace

namespace spz {
namespace python {

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
        .value("SPZ_ADOBE_environment_map", spz::SpzExtensionType::SPZ_ADOBE_environment_map,
               "Adobe environment map extension — an equirectangular RGB map with its frame")
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
        never converted. name is a free-form label. Saving fails unless names are 1-65535 bytes of printable
        ASCII and unique in the cloud, bits is 32 (lossless float32) or 1-16 (quantized per component), values
        holds num_points * dim floats, point-major, and the file version is 4 or later.)doc")
        .def(nb::init<>())
        .def_rw("name", &spz::SpzExtensionCustomAttributeNiantic::name, "Free-form label of printable ASCII, unique per cloud")
        .def_rw("semantic", &spz::SpzExtensionCustomAttributeNiantic::semantic,
                "SPATIAL values are converted with the cloud's coordinates; CUSTOM never are")
        .def_rw("dim", &spz::SpzExtensionCustomAttributeNiantic::dim, "Floats per point (1-65535)")
        .def_rw("bits", &spz::SpzExtensionCustomAttributeNiantic::bits,
                "32 for float32, or 1-16 for quantized storage")
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

    using EnvironmentMap = spz::SpzExtensionEnvironmentMapAdobe;
    using Vec3 = std::array<float, 3>;
    nb::class_<EnvironmentMap, spz::SpzExtensionBase>(m, "SpzExtensionEnvironmentMapAdobe",
        "Equirectangular RGB environment map with its frame; at most one per cloud.")
        .def(nb::init<>())
        .def("__init__",
             [](EnvironmentMap* self, uint32_t width, uint32_t height, const FloatArray& pixels, Vec3 origin,
                Vec3 up, Vec3 e1, Vec3 e2, std::string metadata) {
                 new (self) EnvironmentMap();
                 self->width = width;
                 self->height = height;
                 self->pixels.assign(pixels.data(), pixels.data() + pixels.size());
                 self->origin = origin;
                 self->up = up;
                 self->e1 = e1;
                 self->e2 = e2;
                 self->metadata = std::move(metadata);
             },
             nb::arg("width"), nb::arg("height"), nb::arg("pixels"), nb::arg("origin") = Vec3{0, 0, 0},
             nb::arg("up") = Vec3{0, 1, 0}, nb::arg("e1") = Vec3{1, 0, 0}, nb::arg("e2") = Vec3{0, 0, 1},
             nb::arg("metadata") = "")
        .def_rw("width", &EnvironmentMap::width)
        .def_rw("height", &EnvironmentMap::height)
        .def_prop_rw("pixels",
                     [](const EnvironmentMap& self) { return copyToNumpy(self.pixels); },
                     [](EnvironmentMap& self, const FloatArray& pixels) {
                         self.pixels.assign(pixels.data(), pixels.data() + pixels.size());
                     },
                     nb::rv_policy::move,
                     "Flat float32 array of height * width * 3 values, row-major, RGB interleaved")
        .def_rw("origin", &EnvironmentMap::origin, "Lookup origin, converted like positions")
        .def_rw("up", &EnvironmentMap::up, "Pole direction; pixel row = angle from up")
        .def_rw("e1", &EnvironmentMap::e1, "Azimuth zero direction")
        .def_rw("e2", &EnvironmentMap::e2, "Azimuth +90 degree direction")
        .def_rw("metadata", &EnvironmentMap::metadata, "Free-form UTF-8 (JSON by convention), never converted")
        .def_static("type", &EnvironmentMap::type, "Static method to get the extension type enum value");
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

