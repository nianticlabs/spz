// Extension bindings for the spz Python module.
// This file contains all Python bindings related to SPZ extensions.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>

#include "src/cc/load-spz.h"
#include "src/cc/splat-types.h"
#include "extensions/cc/splat-extensions.h"
#include "extensions/cc/safe-orbit-camera-adobe.h"
#include "extensions/cc/coordinate-system-adobe.h"
#include "extensions/cc/custom-attribute-adobe.h"
#include "extensions/cc/environment-map-adobe.h"
#include "extensions/python/splat-extensions.h"

namespace nb = nanobind;

using FloatArray = nb::ndarray<nb::numpy, float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

namespace spz {
namespace python {
namespace {
FloatArray toNumpy(const std::vector<float>& values) {
    auto* copy = new std::vector<float>(values);
    nb::capsule owner(copy, [](void* p) noexcept { delete static_cast<std::vector<float>*>(p); });
    return FloatArray(copy->data(), {copy->size()}, owner);
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
        .value("SPZ_ADOBE_custom_attribute", spz::SpzExtensionType::SPZ_ADOBE_custom_attribute,
               "Adobe custom attribute extension — one named per-point attribute")
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

    using Attribute = spz::SpzExtensionCustomAttributeAdobe;
    nb::enum_<spz::SpzCustomAttributeKind>(m, "SpzCustomAttributeKind",
                                           "How a custom attribute follows coordinate-system conversions.")
        .value("SCALAR", spz::SpzCustomAttributeKind::Scalar, "Never converted")
        .value("VECTOR3", spz::SpzCustomAttributeKind::Vector3,
               "Each (x, y, z) triple is a direction, converted like positions");
    nb::class_<Attribute, spz::SpzExtensionBase>(m, "SpzExtensionCustomAttributeAdobe")
        .def(nb::init<>())
        .def("__init__",
             [](Attribute* self, std::string name, const FloatArray& values, uint16_t dim,
                spz::SpzCustomAttributeKind kind, uint8_t bits) {
                 new (self) Attribute();
                 self->name = std::move(name);
                 self->values.assign(values.data(), values.data() + values.shape(0));
                 self->dim = dim;
                 self->kind = kind;
                 self->bits = bits;
             },
             nb::arg("name"), nb::arg("values"), nb::arg("dim") = 1,
             nb::arg("kind") = spz::SpzCustomAttributeKind::Scalar, nb::arg("bits") = 32)
        .def_rw("name", &Attribute::name, "Attribute name, unique within a cloud")
        .def_rw("dim", &Attribute::dim, "Values per point")
        .def_rw("kind", &Attribute::kind, "How the values follow coordinate-system conversions")
        .def_rw("bits", &Attribute::bits,
                "32 stores float32 exactly; 1-8 quantizes each component to 2^bits levels")
        .def_prop_rw("values",
                     [](const Attribute& self) { return toNumpy(self.values); },
                     [](Attribute& self, const FloatArray& values) {
                         self.values.assign(values.data(), values.data() + values.shape(0));
                     },
                     nb::rv_policy::move,
                     "Flat float32 array of num_points * dim values, point-major")
        .def_static("type", &Attribute::type, "Static method to get the extension type enum value");

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
                 self->pixels.assign(pixels.data(), pixels.data() + pixels.shape(0));
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
                     [](const EnvironmentMap& self) { return toNumpy(self.pixels); },
                     [](EnvironmentMap& self, const FloatArray& pixels) {
                         self.pixels.assign(pixels.data(), pixels.data() + pixels.shape(0));
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

