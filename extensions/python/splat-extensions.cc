// Extension bindings for the spz Python module.
// This file contains all Python bindings related to SPZ extensions.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>

#include "src/cc/load-spz.h"
#include "src/cc/splat-types.h"
#include "extensions/cc/splat-extensions.h"
#include "extensions/cc/safe-orbit-camera-adobe.h"
#include "extensions/cc/coordinate-system-adobe.h"
#include "extensions/cc/custom-attribute-adobe.h"
#include "extensions/python/splat-extensions.h"

namespace nb = nanobind;

using FloatArray = nb::ndarray<nb::numpy, float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

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
        .value("SPZ_ADOBE_custom_attribute", spz::SpzExtensionType::SPZ_ADOBE_custom_attribute,
               "Adobe custom attribute extension — one named per-point attribute")
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
                     [](const Attribute& self) {
                         auto* copy = new std::vector<float>(self.values);
                         nb::capsule owner(copy, [](void* p) noexcept {
                             delete static_cast<std::vector<float>*>(p);
                         });
                         return FloatArray(copy->data(), {copy->size()}, owner);
                     },
                     [](Attribute& self, const FloatArray& values) {
                         self.values.assign(values.data(), values.data() + values.shape(0));
                     },
                     nb::rv_policy::move,
                     "Flat float32 array of num_points * dim values, point-major")
        .def_static("type", &Attribute::type, "Static method to get the extension type enum value");
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

