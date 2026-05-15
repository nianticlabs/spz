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
#include <cstdlib>
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

// WASM-side unpack chunk size. 64K points keeps the scratch ~18 MB at SH=4.
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

// Build the EmGaussianCloud from already-parsed packed data. `packed` is consumed
// attribute-by-attribute during streaming, so its buffers are released as soon
// as each attribute has been written into a JS-side Float32Array.
EmGaussianCloud buildEmCloud(spz::PackedGaussians packed, const spz::UnpackOptions& options) {
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
  // memory) so they don't count against the WASM heap budget. For very large
  // splats at SH degree 3+, the single SH Float32Array can hit browser
  // TypedArray size limits — those callers should use `loadSpzStreaming` /
  // `loadSpzStreamingAsync` instead, which routes chunks directly into
  // caller-managed destinations.
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

EmGaussianCloud loadSpzFromBuffer(const emscripten::val& buffer, const spz::UnpackOptions& options) {
  // Allocate one input region in the WASM heap and copy directly from the
  // JS-side Uint8Array via HEAPU8.set. Peak memory while parsing is ~2x the
  // input because the caller's JS Uint8Array remains alive alongside the heap
  // copy — to avoid that, use the streaming API and free the JS reference
  // before driving the decode.
  const size_t byteLength = buffer["byteLength"].as<size_t>();
  uint8_t* data = byteLength > 0 ? static_cast<uint8_t*>(std::malloc(byteLength)) : nullptr;
  if (data != nullptr) {
    emscripten::val::module_property("HEAPU8").call<void>(
        "set", buffer, reinterpret_cast<uintptr_t>(data));
  }

  spz::PackedGaussians packed = spz::loadSpzPacked(data, byteLength);
  std::free(data);

  return buildEmCloud(std::move(packed), options);
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
  vectorFromJsArray(emCloud.scales,    cloud.scales);
  vectorFromJsArray(emCloud.rotations, cloud.rotations);
  vectorFromJsArray(emCloud.alphas,    cloud.alphas);
  vectorFromJsArray(emCloud.colors,    cloud.colors);
  vectorFromJsArray(emCloud.sh,        cloud.sh);

  std::vector<uint8_t> output;
  if (!spz::saveSpz(cloud, options, &output)) {
    return emscripten::val::null();
  }

  return jsUint8ArrayFromVector(output);
}

// Streaming API: invoke JS callbacks per attribute chunk so the caller can
// route each chunk directly into its final destination (GPU buffer, multiple
// smaller arrays, OPFS, etc.) without ever allocating a JS-heap buffer that
// scales with numPoints.
//
// Callbacks object:
//   onHeader?(info: { numPoints, shDegree, antialiased, extensions })
//   onChunk(attr: SplatAttribute, pointOffset, data: Float32Array)
//   onDone?()
//   onError?(message: string)
//
// `data` in onChunk is a TypedArray view over WASM linear memory and is valid
// only for the duration of the call. Consumers must copy or upload it before
// returning — a later WASM allocation may detach the view.
void loadSpzStreaming(uintptr_t ptr, size_t byteLength,
                      const spz::UnpackOptions& options,
                      const emscripten::val& callbacks) {
  emscripten::val onHeader = callbacks["onHeader"];
  emscripten::val onChunk  = callbacks["onChunk"];
  emscripten::val onDone   = callbacks["onDone"];
  emscripten::val onError  = callbacks["onError"];

  auto isCallable = [](const emscripten::val& v) {
    return !v.isUndefined() && !v.isNull();
  };
  auto reportError = [&](const std::string& msg) {
    if (isCallable(onError)) onError(msg);
  };

  if (!isCallable(onChunk)) {
    reportError("loadSpzStreaming: callbacks.onChunk is required");
    return;
  }

  const uint8_t* data = reinterpret_cast<const uint8_t*>(ptr);
  spz::PackedGaussians packed = spz::loadSpzPacked(data, byteLength);
  if (packed.numPoints == 0) {
    reportError("loadSpzStreaming: parse returned no points");
    return;
  }

  if (isCallable(onHeader)) {
    emscripten::val header = emscripten::val::object();
    header.set("numPoints", packed.numPoints);
    header.set("shDegree", packed.shDegree);
    header.set("antialiased", packed.antialiased);
#ifdef SPZ_BUILD_EXTENSIONS
    header.set("extensions", jsArrayFromVector(packed.extensions));
#else
    header.set("extensions", emscripten::val::array());
#endif
    onHeader(header);
  }

  const size_t maxFpp = std::max<size_t>(
      4,
      spz::floatsPerPoint(spz::SplatAttribute::Sh, packed.shDegree));
  std::vector<float> scratch(static_cast<size_t>(kStreamingChunkPoints) * maxFpp);

  for (size_t ai = 0; ai < spz::kAllSplatAttributes.size(); ++ai) {
    const spz::SplatAttribute attr = spz::kAllSplatAttributes[ai];
    const size_t fpp = spz::floatsPerPoint(attr, packed.shDegree);
    if (fpp == 0) {
      spz::releasePacked(packed, attr);
      continue;
    }
    for (int32_t off = 0; off < packed.numPoints; off += kStreamingChunkPoints) {
      const int32_t cnt = std::min(kStreamingChunkPoints, packed.numPoints - off);
      const size_t floats = static_cast<size_t>(cnt) * fpp;
      if (!spz::unpackChunk(packed, attr, off, cnt, scratch.data(), options.to)) {
        spz::releasePacked(packed, attr);
        reportError("loadSpzStreaming: unpackChunk failed");
        return;
      }
      // View over WASM scratch — zero-copy. Detaches on any subsequent WASM
      // allocation that grows the heap; the JS callback must consume now.
      onChunk(attr, off,
              emscripten::typed_memory_view(floats, scratch.data()));
    }
    spz::releasePacked(packed, attr);
  }

  if (isCallable(onDone)) onDone();
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

  // Attribute tag passed to loadSpzStreaming's onChunk callback. Registered as
  // `number` so JS receives plain integers (0..5), matching the TS enum in the
  // generated .d.ts. With the default `object` policy, JS would receive enum
  // wrapper instances, and `attr === SplatAttribute.Sh` numeric comparison
  // would fail.
  emscripten::enum_<spz::SplatAttribute>("SplatAttribute",
                                         emscripten::enum_value_type::number)
      .value("Positions", spz::SplatAttribute::Positions)
      .value("Alphas",    spz::SplatAttribute::Alphas)
      .value("Colors",    spz::SplatAttribute::Colors)
      .value("Scales",    spz::SplatAttribute::Scales)
      .value("Rotations", spz::SplatAttribute::Rotations)
      .value("Sh",        spz::SplatAttribute::Sh);

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
  emscripten::function("loadSpzStreaming", &loadSpzStreaming);
  emscripten::function("saveSpzToBuffer", &saveSpzToBuffer);
  emscripten::function("SpzHasExtensionSupport", &spz::hasExtensionSupport);
  emscripten::constant("LATEST_SPZ_HEADER_VERSION", spz::LATEST_SPZ_HEADER_VERSION);
  emscripten::constant("STREAM_CHUNK_POINTS", kStreamingChunkPoints);
}
