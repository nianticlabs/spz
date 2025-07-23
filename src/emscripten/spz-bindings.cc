#include <emscripten/bind.h>

#include <array>
#include <cstdint>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "load-spz.h"

// This file provides Emscripten bindings for the SPZ library, allowing
// JavaScript to interact with SPZ Gaussian splats. A design rule is that
// we avoid exposing C++ STL vectors directly to JavaScript, instead we
// wrap all the functions by using JavaScript-friendly types so that we
// keep the C++ memory managed by Emscripten and avoid memory leaks.

namespace {
// A counterpart to the C++ GaussianCloud structure, used to
// transfer data between C++ and JavaScript.
struct EmGaussianCloud {
  int32_t         numPoints;
  int32_t         shDegree;
  bool            antialiased;
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

template <typename T>
inline void vectorFromJsArray(const emscripten::val& array, std::vector<T>& out) {
  const size_t length = array["length"].as<size_t>();
  out.resize(length);
  for (size_t i = 0; i < length; ++i) {
    out[i] = array[i].as<T>();
  }
}

EmGaussianCloud loadSpzFromBuffer(const emscripten::val& buffer, const spz::UnpackOptions& options) {
  std::vector<uint8_t> bufferInternal;
  vectorFromJsArray(buffer, bufferInternal);

  spz::GaussianCloud cloud = spz::loadSpz(bufferInternal, options);

  EmGaussianCloud emCloud;
  emCloud.numPoints = cloud.numPoints;
  emCloud.shDegree = cloud.shDegree;
  emCloud.antialiased = cloud.antialiased;
  emCloud.positions = jsFloat32ArrayFromVector(cloud.positions);
  emCloud.scales = jsFloat32ArrayFromVector(cloud.scales);
  emCloud.rotations = jsFloat32ArrayFromVector(cloud.rotations);
  emCloud.alphas = jsFloat32ArrayFromVector(cloud.alphas);
  emCloud.colors = jsFloat32ArrayFromVector(cloud.colors);
  emCloud.sh = jsFloat32ArrayFromVector(cloud.sh);

  return emCloud;
}

emscripten::val saveSpzToBuffer(const EmGaussianCloud& emCloud, const spz::PackOptions& options) {
  spz::GaussianCloud cloud;
  cloud.numPoints = emCloud.numPoints;
  cloud.shDegree = emCloud.shDegree;
  cloud.antialiased = emCloud.antialiased;

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
      .value("RUF", spz::CoordinateSystem::RUF);

  emscripten::value_object<spz::PackOptions>("PackOptions")
      .field("from", &spz::PackOptions::from)
      .field("sh1Bits", &spz::PackOptions::sh1Bits)
      .field("shRestBits", &spz::PackOptions::shRestBits)
      .field("hasSafeOrbit", &spz::PackOptions::hasSafeOrbit)
      .field("safeOrbitElevationMin", &spz::PackOptions::safeOrbitElevationMin)
      .field("safeOrbitElevationMax", &spz::PackOptions::safeOrbitElevationMax)
      .field("safeOrbitRadiusMin", &spz::PackOptions::safeOrbitRadiusMin);

  emscripten::value_object<spz::UnpackOptions>("UnpackOptions").field("to", &spz::UnpackOptions::to);

  emscripten::value_object<EmGaussianCloud>("GaussianCloud")
      .field("numPoints", &EmGaussianCloud::numPoints)
      .field("shDegree", &EmGaussianCloud::shDegree)
      .field("antialiased", &EmGaussianCloud::antialiased)
      .field("positions", &EmGaussianCloud::positions)
      .field("scales", &EmGaussianCloud::scales)
      .field("rotations", &EmGaussianCloud::rotations)
      .field("alphas", &EmGaussianCloud::alphas)
      .field("colors", &EmGaussianCloud::colors)
      .field("sh", &EmGaussianCloud::sh);

  emscripten::function("loadSpzFromBuffer", &loadSpzFromBuffer);
  emscripten::function("saveSpzToBuffer", &saveSpzToBuffer);
}
