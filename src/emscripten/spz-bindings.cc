/*
MIT License

Copyright (c) 2025 Adobe Inc.

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

#include <emscripten/bind.h>

#include <array>
#include <cstdint>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "load-spz.h"
#include "utils.h"
#ifdef SPZ_BUILD_EXTENSIONS
#include "extensions/emscripten/splat-extensions.h"
#endif

// This file provides Emscripten bindings for the SPZ library, allowing
// JavaScript to interact with SPZ Gaussian splats. A design rule is that
// we avoid exposing C++ STL vectors directly to JavaScript, instead we
// wrap all the functions by using JavaScript-friendly types so that we
// keep the C++ memory managed by Emscripten and avoid memory leaks.

namespace {
// A counterpart to the C++ GaussianCloud structure, used to
// transfer data between C++ and JavaScript.
struct EmGaussianCloud {
  int32_t numPoints;
  int32_t shDegree;
  bool antialiased;
  emscripten::val extensions;  // Always present, but empty array when extensions disabled
  emscripten::val positions;
  emscripten::val scales;
  emscripten::val rotations;
  emscripten::val alphas;
  emscripten::val colors;
  emscripten::val sh;
};

// A counterpart to the C++ PackOptions structure, used to
// transfer data between C++ and JavaScript.
struct EmPackOptions {
  uint32_t version;
  spz::CoordinateSystem from;
  uint8_t sh1Bits;
  uint8_t shRestBits;
};

inline emscripten::val jsUint8ArrayFromVector(const std::vector<uint8_t>& buffer) {
  emscripten::val Uint8Array = emscripten::val::global("Uint8Array");
  emscripten::val array = Uint8Array.new_(buffer.size());
  for (size_t i = 0; i < buffer.size(); ++i) {
    array.set(i, buffer[i]);
  }
  return array;
}

inline emscripten::val jsFloat32ArrayFromVector(const std::vector<float>& buffer) {
  emscripten::val Float32Array = emscripten::val::global("Float32Array");
  emscripten::val array = Float32Array.new_(buffer.size());
  for (size_t i = 0; i < buffer.size(); ++i) {
    array.set(i, buffer[i]);
  }
  return array;
}

EmGaussianCloud loadSpzFromBuffer(const emscripten::val& buffer, const spz::UnpackOptions& options) {
  std::vector<uint8_t> bufferInternal;
  vectorFromJsArray(buffer, bufferInternal);

  spz::GaussianCloud cloud = spz::loadSpz(bufferInternal, options);

  EmGaussianCloud emCloud;
  emCloud.numPoints = cloud.numPoints;
  emCloud.shDegree = cloud.shDegree;
  emCloud.antialiased = cloud.antialiased;
#ifdef SPZ_BUILD_EXTENSIONS
  emCloud.extensions = jsArrayFromVector(cloud.extensions);
#else
  emCloud.extensions = emscripten::val::array();  // Empty array when extensions disabled
#endif
  emCloud.positions = jsFloat32ArrayFromVector(cloud.positions);
  emCloud.scales = jsFloat32ArrayFromVector(cloud.scales);
  emCloud.rotations = jsFloat32ArrayFromVector(cloud.rotations);
  emCloud.alphas = jsFloat32ArrayFromVector(cloud.alphas);
  emCloud.colors = jsFloat32ArrayFromVector(cloud.colors);
  emCloud.sh = jsFloat32ArrayFromVector(cloud.sh);

  return emCloud;
}

emscripten::val saveSpzToBuffer(const EmGaussianCloud& emCloud, const EmPackOptions& emOptions) {
  spz::GaussianCloud cloud;
  cloud.numPoints = emCloud.numPoints;
  cloud.shDegree = emCloud.shDegree;
  cloud.antialiased = emCloud.antialiased;

#ifdef SPZ_BUILD_EXTENSIONS
  vectorFromJsArray(emCloud.extensions, cloud.extensions);
#endif
  vectorFromJsArray(emCloud.positions, cloud.positions);
  vectorFromJsArray(emCloud.scales, cloud.scales);
  vectorFromJsArray(emCloud.rotations, cloud.rotations);
  vectorFromJsArray(emCloud.alphas, cloud.alphas);
  vectorFromJsArray(emCloud.colors, cloud.colors);
  vectorFromJsArray(emCloud.sh, cloud.sh);

  spz::PackOptions options;
  options.version = emOptions.version;
  options.from = emOptions.from;
  options.sh1Bits = emOptions.sh1Bits;
  options.shRestBits = emOptions.shRestBits;

  std::vector<uint8_t> output;
  if (!spz::saveSpz(cloud, options, &output)) {
    return emscripten::val::null();
  }

  return jsUint8ArrayFromVector(output);
}

}  // namespace

EMSCRIPTEN_BINDINGS(spz_module) {
  // Wrap CoordinateSystem enum
  emscripten::enum_<spz::CoordinateSystem>("CoordinateSystem")
      .value("UNSPECIFIED", spz::CoordinateSystem::UNSPECIFIED)
      .value("LDB", spz::CoordinateSystem::LDB)
      .value("RDB", spz::CoordinateSystem::RDB)
      .value("LUB", spz::CoordinateSystem::LUB)
      .value("RUB", spz::CoordinateSystem::RUB)
      .value("LDF", spz::CoordinateSystem::LDF)
      .value("RDF", spz::CoordinateSystem::RDF)
      .value("LUF", spz::CoordinateSystem::LUF)
      .value("RUF", spz::CoordinateSystem::RUF);

#ifdef SPZ_BUILD_EXTENSIONS
  spz::emscripten::register_extensions();
#endif

  emscripten::value_object<EmPackOptions>("PackOptions")
      .field("version", &EmPackOptions::version)
      .field("from", &EmPackOptions::from)
      .field("sh1Bits", &EmPackOptions::sh1Bits)
      .field("shRestBits", &EmPackOptions::shRestBits);

  emscripten::value_object<spz::UnpackOptions>("UnpackOptions").field("to", &spz::UnpackOptions::to);

  emscripten::value_object<EmGaussianCloud>("GaussianCloud")
      .field("numPoints", &EmGaussianCloud::numPoints)
      .field("shDegree", &EmGaussianCloud::shDegree)
      .field("antialiased", &EmGaussianCloud::antialiased)
      .field("extensions", &EmGaussianCloud::extensions)
      .field("positions", &EmGaussianCloud::positions)
      .field("scales", &EmGaussianCloud::scales)
      .field("rotations", &EmGaussianCloud::rotations)
      .field("alphas", &EmGaussianCloud::alphas)
      .field("colors", &EmGaussianCloud::colors)
      .field("sh", &EmGaussianCloud::sh);

  emscripten::function("loadSpzFromBuffer", &loadSpzFromBuffer);
  emscripten::function("saveSpzToBuffer", &saveSpzToBuffer);
  emscripten::function("SpzHasExtensionSupport", &spz::hasExtensionSupport);
  emscripten::constant("LATEST_SPZ_HEADER_VERSION", spz::LATEST_SPZ_HEADER_VERSION);
}
