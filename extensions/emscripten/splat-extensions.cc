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
#include "extensions/cc/custom-attribute-adobe.h"
#include "extensions/emscripten/splat-extensions.h"
#include "src/emscripten/utils.h"

namespace spz {
namespace emscripten {
namespace {
::emscripten::val customAttributeValues(const spz::SpzExtensionCustomAttributeAdobe& attr) {
  return jsFloat32ArrayFromVector(attr.values);
}

void setCustomAttributeValues(spz::SpzExtensionCustomAttributeAdobe& attr, ::emscripten::val values) {
  vectorFromJsArray(values, attr.values);
}
}  // namespace

// Register all extension-related Emscripten bindings
void register_extensions() {
  ::emscripten::enum_<spz::SpzExtensionType>("SpzExtensionType")
      .value("SPZ_ADOBE_safe_orbit_camera", spz::SpzExtensionType::SPZ_ADOBE_safe_orbit_camera)
      .value("SPZ_ADOBE_coordinate_system", spz::SpzExtensionType::SPZ_ADOBE_coordinate_system)
      .value("SPZ_ADOBE_custom_attribute", spz::SpzExtensionType::SPZ_ADOBE_custom_attribute);

  ::emscripten::enum_<spz::SpzCustomAttributeKind>("SpzCustomAttributeKind")
      .value("Scalar", spz::SpzCustomAttributeKind::Scalar)
      .value("Vector3", spz::SpzCustomAttributeKind::Vector3);

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

  ::emscripten::class_<spz::SpzExtensionCustomAttributeAdobe, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionCustomAttributeAdobe")
      .smart_ptr_constructor("SpzExtensionCustomAttributeAdobePtr", &std::make_shared<spz::SpzExtensionCustomAttributeAdobe>)
      .property("name", &spz::SpzExtensionCustomAttributeAdobe::name)
      .property("dim", &spz::SpzExtensionCustomAttributeAdobe::dim)
      .property("kind", &spz::SpzExtensionCustomAttributeAdobe::kind)
      .property("bits", &spz::SpzExtensionCustomAttributeAdobe::bits)
      .property("values", &customAttributeValues, &setCustomAttributeValues)
      .class_function("type", &spz::SpzExtensionCustomAttributeAdobe::type);

  ::emscripten::function("isKnownPlyExtensionElement", &spz::isKnownPlyExtensionElement);
}

}  // namespace emscripten
}  // namespace spz

