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

#include <vector>

#include "src/cc/load-spz.h"
#include "extensions/cc/splat-extensions.h"
#include "extensions/emscripten/splat-extensions.h"
#include "src/emscripten/utils.h"

namespace spz {
namespace emscripten {

// Register all extension-related Emscripten bindings
void register_extensions() {
  ::emscripten::enum_<spz::SpzExtensionType>("SpzExtensionType")
      .value("SPZ_ADOBE_sh_quantization", spz::SpzExtensionType::SPZ_ADOBE_sh_quantization)
      .value("SPZ_ADOBE_safe_orbit_camera", spz::SpzExtensionType::SPZ_ADOBE_safe_orbit_camera);

  ::emscripten::class_<spz::SpzExtensionBase>("SpzExtensionBase")
      .smart_ptr<std::shared_ptr<spz::SpzExtensionBase>>("SpzExtensionBasePtr")
      .property("extensionType", &spz::SpzExtensionBase::extensionType);

  ::emscripten::class_<spz::SpzExtensionSHQuantizationAdobe, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionSHQuantizationAdobe")
      .constructor<>()
      .property("sh1Bits", &spz::SpzExtensionSHQuantizationAdobe::sh1Bits)
      .property("shRestBits", &spz::SpzExtensionSHQuantizationAdobe::shRestBits)
      .class_function("type", &spz::SpzExtensionSHQuantizationAdobe::type);

  ::emscripten::class_<spz::SpzExtensionSafeOrbitCameraAdobe, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionSafeOrbitCameraAdobe")
      .constructor<>()
      .property("safeOrbitElevationMin", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMin)
      .property("safeOrbitElevationMax", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMax)
      .property("safeOrbitRadiusMin", &spz::SpzExtensionSafeOrbitCameraAdobe::safeOrbitRadiusMin)
      .class_function("type", &spz::SpzExtensionSafeOrbitCameraAdobe::type);
}

}  // namespace emscripten
}  // namespace spz

