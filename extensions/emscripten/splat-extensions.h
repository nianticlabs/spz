// Extension bindings header for the spz Emscripten module.

#ifndef SPLAT_EMSCRIPTEN_EXTENSIONS_H_
#define SPLAT_EMSCRIPTEN_EXTENSIONS_H_

#include <emscripten/val.h>
#include <emscripten/bind.h>

// Forward declarations
namespace spz {
    struct GaussianCloud;
    struct PackOptions;
}

namespace spz {
namespace emscripten {

// Helper function to set extensions in EmGaussianCloud from GaussianCloud
::emscripten::val setExtensionsFromCloud(const spz::GaussianCloud& cloud);

// Helper function to set extensions in GaussianCloud from EmGaussianCloud
void setExtensionsToCloud(const ::emscripten::val& emExtensions, spz::GaussianCloud& cloud);

// Helper function to get extensions from PackOptions
::emscripten::val getExtensionsFromPackOptions(const spz::PackOptions& options);

// Helper function to set extensions in PackOptions
void setExtensionsToPackOptions(const ::emscripten::val& emExtensions, spz::PackOptions& options);

// Register all extension-related Emscripten bindings
void register_extensions();

// Register extension-related fields for PackOptions
void register_pack_options_extensions(::emscripten::value_object<spz::PackOptions>& pack_options);

}  // namespace emscripten
}  // namespace spz

#endif  // SPLAT_EMSCRIPTEN_EXTENSIONS_H_

