// Extension bindings for the spz Python module.
// This file contains all Python bindings related to SPZ extensions.

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>

#include "src/cc/load-spz.h"
#include "src/cc/splat-types.h"
#include "extensions/cc/splat-extensions.h"
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
        .value("SPZ_ADOBE_sh_quantization", spz::SpzExtensionType::SPZ_ADOBE_sh_quantization,
               "Adobe spherical harmonics quantization extension")
        .value("SPZ_ADOBE_safe_orbit_camera", spz::SpzExtensionType::SPZ_ADOBE_safe_orbit_camera,
               "Adobe safe orbit camera extension")
        .export_values();

    nb::class_<spz::SpzExtensionBase>(m, "SpzExtensionBase")
        .def_ro("extension_type", &spz::SpzExtensionBase::extensionType,
                "Type of the SPZ extension");

    nb::class_<spz::SpzExtensionSHQuantizationAdobe, spz::SpzExtensionBase>(m, "SpzExtensionSHQuantizationAdobe")
        .def(nb::init<>())
        .def_rw("sh1_bits", &spz::SpzExtensionSHQuantizationAdobe::sh1Bits,
                "Bits used for first-order spherical harmonics quantization")
        .def_rw("sh_rest_bits", &spz::SpzExtensionSHQuantizationAdobe::shRestBits,
                "Bits used for non-first-order spherical harmonics quantization")
        .def_static("type", &spz::SpzExtensionSHQuantizationAdobe::type,
                    "Static method to get the extension type enum value");

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
}

// Register extension-related properties for PackOptions
void register_pack_options_extensions(nb::class_<spz::PackOptions>& pack_options) {
    pack_options
        .def_rw("extensions", &spz::PackOptions::extensions,
                "List of SPZ extensions associated with this PackOptions. "
                "Users should construct SpzExtensionSHQuantizationAdobe and push it into this vector to enable SH quantization.");
}

// Register extension-related properties for GaussianCloud
void register_gaussian_cloud_extensions(nb::class_<spz::GaussianCloud>& cloud) {
    cloud.def_rw("extensions", &spz::GaussianCloud::extensions,
                 "List of SPZ extensions associated with this GaussianCloud.");
}

}  // namespace python
}  // namespace spz

