/*
MIT License

Copyright (c) 2026 Adobe Inc.

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

#ifndef SPZ_ENVIRONMENT_MAP_ADOBE_H_
#define SPZ_ENVIRONMENT_MAP_ADOBE_H_

#include <array>

#include "splat-extensions.h"

namespace spz {

// An equirectangular RGB environment map with the frame it is looked up in, e.g. for reflections.
// Pixel column = azimuth about `up` measured from `e1` toward `e2`, -pi..pi -> 0..width; pixel row =
// angle from `up`, 0..pi -> 0..height. `origin`, `up`, `e1` and `e2` follow coordinate conversions
// like positions, so lookups stay correct in any coordinate system. `metadata` is free-form UTF-8
// (JSON by convention) and never converted. At most one per cloud: saving fails with two, and
// loading drops every copy when a file has more than one.
//
// `pixels` holds height * width * 3 floats, row-major, RGB interleaved. In SPZ (version 4 only) they
// are stored as float16 in extension frames after the core streams, so float16 inputs round-trip
// exactly and the record in the extension zone holds only the size, frame and metadata; PLY stores
// float32.
struct SpzExtensionEnvironmentMapAdobe : public SpzExtensionBase,
                                         public SpzFramedExtension,
                                         public SpzKeyedExtension,
                                         public SpzCoordinateDependentExtension {
  uint32_t width = 0;
  uint32_t height = 0;
  std::array<float, 3> origin = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> up = {0.0f, 1.0f, 0.0f};
  std::array<float, 3> e1 = {1.0f, 0.0f, 0.0f};
  std::array<float, 3> e2 = {0.0f, 0.0f, 1.0f};
  std::string metadata;
  std::vector<float> pixels;

  static const std::unordered_set<std::string> kRequiredPlyElementNames;

  SpzExtensionEnvironmentMapAdobe();
  // The pixels only fit in extension frames, which saveSpz writes through planFrameWrites(); write(os)
  // therefore writes nothing (and logs an error), and payloadBytes() is 0.
  uint32_t payloadBytes() const override;
  void write(std::ostream& os) const override;
  SpzExtensionBase* copyAsRawData() const override;
  std::optional<std::shared_ptr<SpzExtensionBase>> tryReadFromPly(
      std::istream& in, const std::vector<PlyExtraElement>& elements) const override;
  void writePlyHeader(std::ostream& out) const override;
  void writePlyData(std::ostream& out) const override;
  static SpzExtensionType type();

  // False (after logging why) if this map can't be saved.
  bool isValid() const override;
  // Every environment map shares one key, so a cloud holds at most one.
  std::string uniquenessKey() const override;
  // A copy with origin, up, e1 and e2 converted like positions, or nullptr if c leaves them unchanged.
  SpzExtensionBasePtr withConvertedCoordinates(const CoordinateConverter& c) const override;
  std::optional<SpzFrameWritePlan> planFrameWrites(int32_t numPoints) const override;
};

}  // namespace spz

#endif  // SPZ_ENVIRONMENT_MAP_ADOBE_H_
