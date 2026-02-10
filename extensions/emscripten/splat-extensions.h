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
// Register all extension-related Emscripten bindings
void register_extensions();
}  // namespace emscripten
}  // namespace spz

#endif  // SPLAT_EMSCRIPTEN_EXTENSIONS_H_

