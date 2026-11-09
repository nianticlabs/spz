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
#include "extensions/cc/temporal-niantic.h"
#include "extensions/emscripten/splat-extensions.h"
#include "src/emscripten/utils.h"

namespace spz {
namespace emscripten {
namespace {
using ::emscripten::val;

// Returns the bytes of a Uint8Array or ArrayBuffer.
std::vector<uint8_t> bytesFromJs(const val& buffer) {
  const val bytes = buffer.instanceof(val::global("ArrayBuffer"))
                        ? val::global("Uint8Array").new_(buffer)
                        : buffer;
  std::vector<uint8_t> out(bytes["byteLength"].as<size_t>());
  val(::emscripten::typed_memory_view(out.size(), out.data())).call<void>("set", bytes);
  return out;
}

val channelToJs(const SpzTemporalChannel& c) {
  val o = val::object();
  o.set("target", static_cast<int>(c.target));
  o.set("type", static_cast<int>(c.type));
  o.set("flags", static_cast<int>(c.flags));
  o.set("dim", c.dim);
  o.set("encoding", static_cast<int>(c.encoding));
  o.set("bits", c.bits);
  o.set("steps", jsFloat32ArrayFromVector(c.steps));
  if (c.pointRange) {
    val range = val::array();
    range.call<void>("push", c.pointRange->begin);
    range.call<void>("push", c.pointRange->end);
    o.set("pointRange", range);
  } else {
    o.set("pointRange", val::null());
  }
  o.set("values", jsFloat32ArrayFromVector(c.values));
  return o;
}

SpzTemporalChannel channelFromJs(const val& o) {
  SpzTemporalChannel c;
  c.target = static_cast<SpzTemporalTarget>(o["target"].as<int>());
  c.type = static_cast<SpzTemporalChannelType>(o["type"].as<int>());
  c.flags = o["flags"].isUndefined() ? SpzTemporalChannelFlags::Required
                                     : static_cast<SpzTemporalChannelFlags>(o["flags"].as<int>());
  c.dim = o["dim"].as<uint16_t>();
  if (!o["encoding"].isUndefined()) c.encoding = static_cast<SpzChannelEncoding>(o["encoding"].as<int>());
  if (!o["bits"].isUndefined()) c.bits = o["bits"].as<uint8_t>();
  if (!o["steps"].isUndefined() && !o["steps"].isNull()) vectorFromJsArray(o["steps"], c.steps);
  const val range = o["pointRange"];
  if (!range.isUndefined() && !range.isNull()) {
    c.pointRange = SpzPointRange{range[0].as<uint32_t>(), range[1].as<uint32_t>()};
  }
  vectorFromJsArray(o["values"], c.values);
  return c;
}

// Plays an animated clip: holds the cloud in the WASM heap and evaluates frames into reused
// buffers. The views returned by positions(), rotations() and alphas() alias WASM memory: they
// are valid until the next evaluate() call or until the heap grows, so copy or upload them
// right away.
class TemporalPlayer {
 public:
  TemporalPlayer(const val& buffer, CoordinateSystem to) {
    UnpackOptions options;
    options.to = to;
    options.temporal = TemporalLoad::Include;
    const std::vector<uint8_t> bytes = bytesFromJs(buffer);
    cloud_ = loadSpz(bytes, options);
    if (cloud_.numPoints == 0) {
      val::global("Error").new_(val("TemporalPlayer: not a valid SPZ file")).throw_();
    }
    model_ = bindTemporalModel(cloud_);
    if (model_) {
      index_ = computeTemporalTimeIndex(*model_, cloud_.numPoints, temporalTimeBinBoundaries(*model_));
    }
    evaluate(model_ ? model_->clip->posterTime : 0.0);
  }

  bool isAnimated() const { return model_.has_value(); }
  int32_t numPoints() const { return cloud_.numPoints; }
  int32_t shDegree() const { return cloud_.shDegree; }
  double clipStart() const { return model_ ? model_->clip->clipStart : 0.0; }
  double clipEnd() const { return model_ ? model_->clip->clipEnd : 0.0; }
  double posterTime() const { return model_ ? model_->clip->posterTime : 0.0; }
  int loop() const { return model_ ? static_cast<int>(model_->clip->loop) : 0; }
  float nominalFps() const { return model_ ? model_->clip->nominalFps : 0.0f; }
  float opacityCutoff() const { return model_ ? model_->clip->opacityCutoff : 0.0f; }

  double clipTime(double playbackSeconds) const {
    return model_ ? spz::clipTime(*model_->clip, playbackSeconds) : 0.0;
  }

  // Evaluates clip time t into the frame buffers. A static cloud always shows its core state.
  void evaluate(double t) {
    if (model_) {
      evaluateTemporalInto(cloud_, *model_, t, &frame_);
    } else {
      frame_.positions = cloud_.positions;
      frame_.rotations = cloud_.rotations;
      frame_.alphas = cloud_.alphas;
    }
  }

  // The most Gaussians evaluateVisible can return: size per-frame buffers for this many.
  uint32_t maxVisible() const {
    return model_ ? index_.maxBinPoints : static_cast<uint32_t>(cloud_.numPoints);
  }

  // Evaluates only the Gaussians that may be visible at clip time t (per the time index) and keeps
  // those at least opacityCutoff opaque; returns how many. A static cloud returns all its points.
  uint32_t evaluateVisible(double t) {
    if (model_) {
      evaluateVisibleTemporal(cloud_, *model_, index_, t, &visible_);
    } else {
      visible_.indices.resize(static_cast<size_t>(cloud_.numPoints));
      for (uint32_t i = 0; i < visible_.indices.size(); ++i) visible_.indices[i] = i;
      visible_.positions = cloud_.positions;
      visible_.rotations = cloud_.rotations;
      visible_.alphas = cloud_.alphas;
    }
    return static_cast<uint32_t>(visible_.indices.size());
  }

  val visibleIndices() const {
    return val(::emscripten::typed_memory_view(visible_.indices.size(), visible_.indices.data()));
  }
  val visiblePositions() const { return view(visible_.positions); }
  val visibleRotations() const { return view(visible_.rotations); }
  val visibleAlphas() const { return view(visible_.alphas); }

  val positions() const { return view(frame_.positions); }
  val rotations() const { return view(frame_.rotations); }
  val alphas() const { return view(frame_.alphas); }
  val scales() const { return view(cloud_.scales); }
  val colors() const { return view(cloud_.colors); }
  val sh() const { return view(cloud_.sh); }

 private:
  static val view(const std::vector<float>& v) {
    return val(::emscripten::typed_memory_view(v.size(), v.data()));
  }

  GaussianCloud cloud_;
  std::optional<SpzTemporalModel> model_;
  SpzTemporalTimeIndex index_;
  SpzTemporalFrame frame_;
  SpzTemporalVisibleFrame visible_;
};

void registerTemporal() {
  using Temporal = spz::SpzExtensionTemporalNiantic;
  // Registered as numbers so JS compares plain integers, matching the TS enums.
  ::emscripten::enum_<spz::TemporalLoad>("TemporalLoad", ::emscripten::enum_value_type::number)
      .value("SKIP", spz::TemporalLoad::Skip)
      .value("INCLUDE", spz::TemporalLoad::Include);
  ::emscripten::enum_<spz::SpzTemporalTarget>("SpzTemporalTarget", ::emscripten::enum_value_type::number)
      .value("POSITION", spz::SpzTemporalTarget::Position)
      .value("ROTATION", spz::SpzTemporalTarget::Rotation)
      .value("OPACITY", spz::SpzTemporalTarget::Opacity)
      .value("TRANSFORM", spz::SpzTemporalTarget::Transform);
  ::emscripten::enum_<spz::SpzTemporalChannelType>("SpzTemporalChannelType", ::emscripten::enum_value_type::number)
      .value("POLYNOMIAL", spz::SpzTemporalChannelType::Polynomial)
      .value("WINDOW", spz::SpzTemporalChannelType::Window)
      .value("SKIN", spz::SpzTemporalChannelType::Skin);
  ::emscripten::enum_<spz::SpzTemporalChannelFlags>("SpzTemporalChannelFlags", ::emscripten::enum_value_type::number)
      .value("NONE", spz::SpzTemporalChannelFlags::None)
      .value("REQUIRED", spz::SpzTemporalChannelFlags::Required);
  ::emscripten::enum_<spz::SpzChannelEncoding>("SpzChannelEncoding", ::emscripten::enum_value_type::number)
      .value("FLOAT32", spz::SpzChannelEncoding::Float32)
      .value("RANGE", spz::SpzChannelEncoding::Range)
      .value("STEP", spz::SpzChannelEncoding::Step);
  ::emscripten::enum_<spz::SpzTemporalLoop>("SpzTemporalLoop", ::emscripten::enum_value_type::number)
      .value("ONCE", spz::SpzTemporalLoop::Once)
      .value("LOOP", spz::SpzTemporalLoop::Loop)
      .value("PING_PONG", spz::SpzTemporalLoop::PingPong);

  ::emscripten::class_<Temporal, ::emscripten::base<spz::SpzExtensionBase>>("SpzExtensionTemporalNiantic")
      .smart_ptr_constructor("SpzExtensionTemporalNiantic", &std::make_shared<Temporal>)
      .property("clipStart", &Temporal::clipStart)
      .property("clipEnd", &Temporal::clipEnd)
      .property("posterTime", &Temporal::posterTime)
      .property("loop", &Temporal::loop)
      .property("nominalFps", &Temporal::nominalFps)
      .property("opacityCutoff", &Temporal::opacityCutoff)
      .property("driverCount", &Temporal::driverCount)
      .function("getChannels", ::emscripten::optional_override([](const Temporal& self) {
        val out = val::array();
        for (const auto& c : self.channels) out.call<void>("push", channelToJs(c));
        return out;
      }))
      .function("setChannels", ::emscripten::optional_override([](Temporal& self, const val& channels) {
        self.channels.clear();
        const size_t count = channels["length"].as<size_t>();
        for (size_t i = 0; i < count; ++i) self.channels.push_back(channelFromJs(channels[i]));
      }))
      .function("getDriverKeyframeTimes", ::emscripten::optional_override([](const Temporal& self) {
        val array = val::global("Float64Array").new_(self.driverKeyframeTimes.size());
        array.call<void>("set", ::emscripten::typed_memory_view(self.driverKeyframeTimes.size(),
                                                                self.driverKeyframeTimes.data()));
        return array;
      }))
      .function("setDriverKeyframeTimes", ::emscripten::optional_override([](Temporal& self, const val& times) {
        vectorFromJsArray(times, self.driverKeyframeTimes);
      }))
      .function("getDriverNames", ::emscripten::optional_override([](const Temporal& self) {
        return jsArrayFromVector(self.driverNames);
      }))
      .function("setDriverNames", ::emscripten::optional_override([](Temporal& self, const val& names) {
        vectorFromJsArray(names, self.driverNames);
      }))
      .function("getDriverTransforms", ::emscripten::optional_override([](const Temporal& self) {
        return jsFloat32ArrayFromVector(self.driverTransforms);
      }))
      .function("setDriverTransforms", ::emscripten::optional_override([](Temporal& self, const val& values) {
        vectorFromJsArray(values, self.driverTransforms);
      }))
      .class_function("type", &Temporal::type);

  ::emscripten::function("clipTime", &spz::clipTime);

  ::emscripten::class_<TemporalPlayer>("TemporalPlayer")
      .constructor<const val&, CoordinateSystem>()
      .function("isAnimated", &TemporalPlayer::isAnimated)
      .function("numPoints", &TemporalPlayer::numPoints)
      .function("shDegree", &TemporalPlayer::shDegree)
      .function("clipStart", &TemporalPlayer::clipStart)
      .function("clipEnd", &TemporalPlayer::clipEnd)
      .function("posterTime", &TemporalPlayer::posterTime)
      .function("loop", &TemporalPlayer::loop)
      .function("nominalFps", &TemporalPlayer::nominalFps)
      .function("opacityCutoff", &TemporalPlayer::opacityCutoff)
      .function("clipTime", &TemporalPlayer::clipTime)
      .function("evaluate", &TemporalPlayer::evaluate)
      .function("maxVisible", &TemporalPlayer::maxVisible)
      .function("evaluateVisible", &TemporalPlayer::evaluateVisible)
      .function("visibleIndices", &TemporalPlayer::visibleIndices)
      .function("visiblePositions", &TemporalPlayer::visiblePositions)
      .function("visibleRotations", &TemporalPlayer::visibleRotations)
      .function("visibleAlphas", &TemporalPlayer::visibleAlphas)
      .function("positions", &TemporalPlayer::positions)
      .function("rotations", &TemporalPlayer::rotations)
      .function("alphas", &TemporalPlayer::alphas)
      .function("scales", &TemporalPlayer::scales)
      .function("colors", &TemporalPlayer::colors)
      .function("sh", &TemporalPlayer::sh);
}
}  // namespace

// Register all extension-related Emscripten bindings
void register_extensions() {
  ::emscripten::enum_<spz::SpzExtensionType>("SpzExtensionType")
      .value("SPZ_ADOBE_safe_orbit_camera", spz::SpzExtensionType::SPZ_ADOBE_safe_orbit_camera)
      .value("SPZ_ADOBE_coordinate_system", spz::SpzExtensionType::SPZ_ADOBE_coordinate_system)
      .value("SPZ_NIANTIC_georeference", spz::SpzExtensionType::SPZ_NIANTIC_georeference)
      .value("SPZ_NIANTIC_custom_attribute", spz::SpzExtensionType::SPZ_NIANTIC_custom_attribute)
      .value("SPZ_NIANTIC_temporal", spz::SpzExtensionType::SPZ_NIANTIC_temporal);

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
      .function("getPointRange", ::emscripten::optional_override([](const CustomAttribute& self) {
        if (!self.pointRange) return val::null();
        val range = val::array();
        range.call<void>("push", self.pointRange->begin);
        range.call<void>("push", self.pointRange->end);
        return range;
      }))
      .function("setPointRange", ::emscripten::optional_override([](CustomAttribute& self, const val& range) {
        if (range.isNull() || range.isUndefined()) {
          self.pointRange.reset();
          return;
        }
        self.pointRange = spz::SpzPointRange{range[0].as<uint32_t>(), range[1].as<uint32_t>()};
      }))
      .class_function("type", &CustomAttribute::type);

  registerTemporal();

  ::emscripten::function("isKnownPlyExtensionElement", &spz::isKnownPlyExtensionElement);
}

}  // namespace emscripten
}  // namespace spz

