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

#include "chunked-unpack.h"
#include "load-spz.h"
#include "splat-utils.h"
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


inline emscripten::val jsUint8ArrayFromVector(const std::vector<uint8_t>& buffer) {
  emscripten::val Uint8Array = emscripten::val::global("Uint8Array");
  emscripten::val array = Uint8Array.new_(buffer.size());
  array.call<void>("set", emscripten::typed_memory_view(buffer.size(), buffer.data()));
  return array;
}

inline emscripten::val jsFloat32ArrayFromVector(const std::vector<float>& buffer) {
  emscripten::val Float32Array = emscripten::val::global("Float32Array");
  emscripten::val array = Float32Array.new_(buffer.size());
  array.call<void>("set", emscripten::typed_memory_view(buffer.size(), buffer.data()));
  return array;
}

// Chunk size of 64K points keeps the scratch buffer ~18 MB at worst.
constexpr int32_t kStreamingChunkPoints = 1 << 16;

// Fill a pre-allocated JS Float32Array with one attribute by chunked unpacking
// through a reusable WASM-side scratch buffer.
inline void streamAttributeIntoJsArray(spz::PackedGaussians& packed,
                                       spz::SplatAttribute attr,
                                       spz::CoordinateSystem coordTo,
                                       std::vector<float>& scratch,
                                       emscripten::val& outArray) {
  const int32_t numPoints = packed.numPoints;
  const size_t fpp = spz::floatsPerPoint(attr, packed.shDegree);
  if (fpp == 0 || numPoints == 0) {
    spz::releasePacked(packed, attr);
    return;
  }
  for (int32_t off = 0; off < numPoints; off += kStreamingChunkPoints) {
    const int32_t cnt = std::min(kStreamingChunkPoints, numPoints - off);
    const size_t floats = static_cast<size_t>(cnt) * fpp;
    if (!spz::unpackChunk(packed, attr, off, cnt, scratch.data(), coordTo)) {
      break;
    }
    outArray.call<void>("set",
                        emscripten::typed_memory_view(floats, scratch.data()),
                        static_cast<size_t>(off) * fpp);
  }
  spz::releasePacked(packed, attr);
}

EmGaussianCloud loadSpzFromBuffer(const emscripten::val& buffer, const spz::UnpackOptions& options) {
  std::vector<uint8_t> bufferInternal;
  vectorFromJsArray(buffer, bufferInternal);

  // NOTE: This can surpass the 2GB WASM heap limit if the input buffer is too large.
  // Consider using a streaming approach if the input buffer is too large.
  spz::PackedGaussians packed = spz::loadSpzPacked(bufferInternal.data(),
                                                   bufferInternal.size());

  // Release the input bytes from the WASM heap.
  std::vector<uint8_t>().swap(bufferInternal);

  EmGaussianCloud emCloud;
  emCloud.numPoints = packed.numPoints;
  emCloud.shDegree = packed.shDegree;
  emCloud.antialiased = packed.antialiased;
#ifdef SPZ_BUILD_EXTENSIONS
  emCloud.extensions = jsArrayFromVector(packed.extensions);
#else
  emCloud.extensions = emscripten::val::array();
#endif

  // Pre-allocate the six output Float32Arrays in JS heap (outside WASM linear
  // memory) so they don't count against the WASM heap budget.
  emscripten::val Float32Array = emscripten::val::global("Float32Array");
  auto allocFloatArray = [&](spz::SplatAttribute attr) {
    return Float32Array.new_(spz::floatCount(packed, attr));
  };
  emCloud.positions = allocFloatArray(spz::SplatAttribute::Positions);
  emCloud.scales    = allocFloatArray(spz::SplatAttribute::Scales);
  emCloud.rotations = allocFloatArray(spz::SplatAttribute::Rotations);
  emCloud.alphas    = allocFloatArray(spz::SplatAttribute::Alphas);
  emCloud.colors    = allocFloatArray(spz::SplatAttribute::Colors);
  emCloud.sh        = allocFloatArray(spz::SplatAttribute::Sh);

  // Reusable WASM-side chunk buffer sized for the widest attribute (SH;
  // floor of 4 covers rotations when shDegree == 0).
  const size_t maxFpp = std::max<size_t>(
      4,
      spz::floatsPerPoint(spz::SplatAttribute::Sh, packed.shDegree));
  std::vector<float> scratch(static_cast<size_t>(kStreamingChunkPoints) * maxFpp);

  emscripten::val* slots[] = {
    &emCloud.positions, &emCloud.alphas, &emCloud.colors,
    &emCloud.scales,    &emCloud.rotations, &emCloud.sh,
  };
  static_assert(sizeof(slots) / sizeof(slots[0]) == spz::kAllSplatAttributes.size(),
                "slot count must match kAllSplatAttributes");
  for (size_t i = 0; i < spz::kAllSplatAttributes.size(); ++i) {
    streamAttributeIntoJsArray(packed, spz::kAllSplatAttributes[i], options.to,
                               scratch, *slots[i]);
  }

  return emCloud;
}

emscripten::val saveSpzToBuffer(const EmGaussianCloud& emCloud, const spz::PackOptions& options) {
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
      .value("RUF", spz::CoordinateSystem::RUF)
      .value("LFD", spz::CoordinateSystem::LFD)
      .value("RFD", spz::CoordinateSystem::RFD)
      .value("LFU", spz::CoordinateSystem::LFU)
      .value("RFU", spz::CoordinateSystem::RFU)
      .value("LBD", spz::CoordinateSystem::LBD)
      .value("RBD", spz::CoordinateSystem::RBD)
      .value("LBU", spz::CoordinateSystem::LBU)
      .value("RBU", spz::CoordinateSystem::RBU);

#ifdef SPZ_BUILD_EXTENSIONS
  spz::emscripten::register_extensions();
#endif

  emscripten::value_object<spz::PackOptions>("PackOptions")
      .field("version", &spz::PackOptions::version)
      .field("from", &spz::PackOptions::from)
      .field("sh1Bits", &spz::PackOptions::sh1Bits)
      .field("shRestBits", &spz::PackOptions::shRestBits);

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
