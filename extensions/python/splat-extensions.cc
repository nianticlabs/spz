// Extension bindings for the spz Python module.
// This file contains all Python bindings related to SPZ extensions.

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>

#include "src/cc/load-spz.h"
#include "src/cc/splat-types.h"
#include "extensions/cc/splat-extensions.h"
#include "extensions/cc/safe-orbit-camera-adobe.h"
#include "extensions/cc/coordinate-system-adobe.h"
#include "extensions/python/splat-extensions.h"

namespace nb = nanobind;

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

