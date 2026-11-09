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
#include "extensions/cc/custom-attribute-niantic.h"
#include "extensions/cc/environment-map-adobe.h"
#include "extensions/emscripten/splat-extensions.h"
#include "src/emscripten/utils.h"

namespace spz {
namespace emscripten {
namespace {
::emscripten::val environmentMapPixels(const spz::SpzExtensionEnvironmentMapAdobe& env) {
  return jsFloat32ArrayFromVector(env.pixels);
}

void setEnvironmentMapPixels(spz::SpzExtensionEnvironmentMapAdobe& env, ::emscripten::val pixels) {
  vectorFromJsArray(pixels, env.pixels);
}
}  // namespace

// Register all extension-related Emscripten bindings
void register_extensions() {
  ::emscripten::enum_<spz::SpzExtensionType>("SpzExtensionType")
      .value("SPZ_ADOBE_safe_orbit_camera", spz::SpzExtensionType::SPZ_ADOBE_safe_orbit_camera)
      .value("SPZ_ADOBE_coordinate_system", spz::SpzExtensionType::SPZ_ADOBE_coordinate_system)
      .value("SPZ_ADOBE_environment_map", spz::SpzExtensionType::SPZ_ADOBE_environment_map)
      .value("SPZ_NIANTIC_georeference", spz::SpzExtensionType::SPZ_NIANTIC_georeference)
      .value("SPZ_NIANTIC_custom_attribute", spz::SpzExtensionType::SPZ_NIANTIC_custom_attribute);

  ::emscripten::value_array<std::array<float, 3>>("Vec3f")
      .element(::emscripten::index<0>())
      .element(::emscripten::index<1>())
      .element(::emscripten::index<2>());

  ::emscripten::class_<spz::SpzExtensionBase>("SpzExtensionBase")
      .smart_ptr<std::shared_ptr<spz::SpzExtensionBase>>("SpzExtensionBasePtr")
      .property("extensionType", &spz::SpzExtensionBase::extensionType);

  // Extensions are held as shared_ptrs (GaussianCloud::extensions); objects made with a plain
  // constructor are raw pointers, which embind refuses to pass to saveSpzToBuffer.
  ::emscripten::class_<spz::SpzExtensionSafeOrbitCameraAdobe, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionSafeOrbitCameraAdobe")
      .smart_ptr_constructor("SpzExtensionSafeOrbitCameraAdobePtr", &std::make_shared<spz::SpzExtensionSafeOrbitCameraAdobe>)
      .property("safeOrbitElevationMin", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMin)
      .property("safeOrbitElevationMax", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMax)
      .property("safeOrbitRadiusMin", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitRadiusMin)
      .class_function("type", &spz::SpzExtensionSafeOrbitCameraAdobe::type);

  ::emscripten::class_<spz::SpzExtensionCoordinateSystemAdobe, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionCoordinateSystemAdobe")
      .smart_ptr_constructor("SpzExtensionCoordinateSystemAdobePtr", &std::make_shared<spz::SpzExtensionCoordinateSystemAdobe>)
      .property("coordinateSystem", &spz::SpzExtensionCoordinateSystemAdobe::coordinateSystem)
      .class_function("type", &spz::SpzExtensionCoordinateSystemAdobe::type);

  // Registered as `number` (as SplatAttribute in spz-bindings.cc) so JS sees plain integers,
  // matching the numeric enum declared in splat-extensions.d.ts.in.
  ::emscripten::enum_<spz::SpzExtensionGeoreferenceNiantic::CrsEncoding>(
      "GeoreferenceCrsEncoding", ::emscripten::enum_value_type::number)
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

  ::emscripten::enum_<spz::SpzAttributeSemantic>("SpzAttributeSemantic")
      .value("CUSTOM", spz::SpzAttributeSemantic::Custom)
      .value("SPATIAL", spz::SpzAttributeSemantic::Spatial);

  using CustomAttribute = spz::SpzExtensionCustomAttributeNiantic;
  // Created through a shared_ptr so attributes built in JS can be added to a cloud's extensions.
  ::emscripten::class_<CustomAttribute, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionCustomAttributeNiantic")
      .smart_ptr_constructor("SpzExtensionCustomAttributeNiantic", &std::make_shared<CustomAttribute>)
      .property("name", &CustomAttribute::name)
      .property("dim", &CustomAttribute::dim)
      .property("semantic", &CustomAttribute::semantic)
      .property("bits", &CustomAttribute::bits)
      .function("getValues", ::emscripten::optional_override([](const CustomAttribute& self) {
        return jsFloat32ArrayFromVector(self.values);
      }))
      .function("setValues", ::emscripten::optional_override([](CustomAttribute& self, const ::emscripten::val& values) {
        vectorFromJsArray(values, self.values);
      }))
      .class_function("type", &CustomAttribute::type);

  using EnvironmentMap = spz::SpzExtensionEnvironmentMapAdobe;
  ::emscripten::class_<EnvironmentMap, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionEnvironmentMapAdobe")
      .smart_ptr_constructor("SpzExtensionEnvironmentMapAdobePtr", &std::make_shared<EnvironmentMap>)
      .property("width", &EnvironmentMap::width)
      .property("height", &EnvironmentMap::height)
      .property("origin", &EnvironmentMap::origin)
      .property("up", &EnvironmentMap::up)
      .property("e1", &EnvironmentMap::e1)
      .property("e2", &EnvironmentMap::e2)
      .property("metadata", &EnvironmentMap::metadata)
      .property("pixels", &environmentMapPixels, &setEnvironmentMapPixels)
      .class_function("type", &EnvironmentMap::type);

  ::emscripten::function("isKnownPlyExtensionElement", &spz::isKnownPlyExtensionElement);
}

}  // namespace emscripten
}  // namespace spz

