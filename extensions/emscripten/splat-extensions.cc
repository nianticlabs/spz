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

#include "src/cc/load-spz.h"
#include "extensions/cc/splat-extensions.h"
#include "extensions/cc/safe-orbit-camera-adobe.h"
#include "extensions/cc/coordinate-system-adobe.h"
#include "extensions/cc/georeference-niantic.h"
#include "extensions/emscripten/splat-extensions.h"
#include "src/emscripten/utils.h"

namespace spz {
namespace emscripten {

// Register all extension-related Emscripten bindings
void register_extensions() {
  ::emscripten::enum_<spz::SpzExtensionType>("SpzExtensionType")
      .value("SPZ_ADOBE_safe_orbit_camera", spz::SpzExtensionType::SPZ_ADOBE_safe_orbit_camera)
      .value("SPZ_ADOBE_coordinate_system", spz::SpzExtensionType::SPZ_ADOBE_coordinate_system)
      .value("SPZ_NIANTIC_georeference", spz::SpzExtensionType::SPZ_NIANTIC_georeference);

  ::emscripten::class_<spz::SpzExtensionBase>("SpzExtensionBase")
      .smart_ptr<std::shared_ptr<spz::SpzExtensionBase>>("SpzExtensionBasePtr")
      .property("extensionType", &spz::SpzExtensionBase::extensionType);

  ::emscripten::class_<spz::SpzExtensionSafeOrbitCameraAdobe, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionSafeOrbitCameraAdobe")
      .constructor<>()
      .property("safeOrbitElevationMin", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMin)
      .property("safeOrbitElevationMax", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMax)
      .property("safeOrbitRadiusMin", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitRadiusMin)
      .class_function("type", &spz::SpzExtensionSafeOrbitCameraAdobe::type);

  ::emscripten::class_<spz::SpzExtensionCoordinateSystemAdobe, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionCoordinateSystemAdobe")
      .constructor<>()
      .property("coordinateSystem", &spz::SpzExtensionCoordinateSystemAdobe::coordinateSystem)
      .class_function("type", &spz::SpzExtensionCoordinateSystemAdobe::type);

  // Registered as `number` (as SplatAttribute in spz-bindings.cc) so JS sees plain integers,
  // matching the numeric enum declared in splat-extensions.d.ts.in.
  ::emscripten::enum_<spz::SpzExtensionGeoreferenceNiantic::CrsEncoding>(
      "CrsEncoding", ::emscripten::enum_value_type::number)
      .value("AUTHORITY_CODE", spz::SpzExtensionGeoreferenceNiantic::CrsEncoding::AuthorityCode)
      .value("PROJJSON", spz::SpzExtensionGeoreferenceNiantic::CrsEncoding::Projjson);

  ::emscripten::value_array<std::array<double, 3>>("ArrayDouble3")
      .element(::emscripten::index<0>())
      .element(::emscripten::index<1>())
      .element(::emscripten::index<2>());

  ::emscripten::value_array<std::array<double, 4>>("ArrayDouble4")
      .element(::emscripten::index<0>())
      .element(::emscripten::index<1>())
      .element(::emscripten::index<2>())
      .element(::emscripten::index<3>());

  ::emscripten::class_<spz::SpzExtensionGeoreferenceNiantic, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionGeoreferenceNiantic")
      .constructor<>()
      .property("crs", &spz::SpzExtensionGeoreferenceNiantic::crs)
      .property("crsEncoding", &spz::SpzExtensionGeoreferenceNiantic::crsEncoding)
      .property("origin", &spz::SpzExtensionGeoreferenceNiantic::origin)
      .property("rotation", &spz::SpzExtensionGeoreferenceNiantic::rotation)
      .property("scale", &spz::SpzExtensionGeoreferenceNiantic::scale)
      .property("epoch", &spz::SpzExtensionGeoreferenceNiantic::epoch)
      .class_function("type", &spz::SpzExtensionGeoreferenceNiantic::type);

  ::emscripten::function("isKnownPlyExtensionElement", &spz::isKnownPlyExtensionElement);
}

}  // namespace emscripten
}  // namespace spz

