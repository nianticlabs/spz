// Extension bindings header for the spz Python module.

#ifndef SPLAT_PYTHON_EXTENSIONS_H_
#define SPLAT_PYTHON_EXTENSIONS_H_

#include <nanobind/nanobind.h>

namespace nb = nanobind;

// Forward declarations
namespace spz {
    class PackOptions;
    class GaussianCloud;
}

namespace spz {
namespace python {

// Register all extension-related Python bindings
void register_extensions(nb::module_& m);

// Register extension-related properties for PackOptions
void register_pack_options_extensions(nb::class_<spz::PackOptions>& pack_options);

// Register extension-related properties for GaussianCloud
void register_gaussian_cloud_extensions(nb::class_<spz::GaussianCloud>& cloud);

}  // namespace python
}  // namespace spz

#endif  // SPLAT_PYTHON_EXTENSIONS_H_

