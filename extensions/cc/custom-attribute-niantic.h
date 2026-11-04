/*
MIT License

Copyright (c) 2026 Niantic Spatial

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

#ifndef SPZ_CUSTOM_ATTRIBUTE_NIANTIC_H_
#define SPZ_CUSTOM_ATTRIBUTE_NIANTIC_H_

#include "splat-extensions.h"

namespace spz {

// How a custom attribute responds to conversions.
enum class SpzAttributeSemantic : uint8_t {
  Custom = 0,   // never converted
  Spatial = 1,  // vectors in the cloud's coordinate space rotated and flipped like positions
};

// A named per-point attribute outside the core set (normals, material parameters, embeddings,
// ...). A cloud may carry any number of them, one record each. The values are stored as zstd
// extension frames after the last core stream, so files stay loadable by every SPZ v4 reader;
// readers without this extension load the core cloud and skip the attribute.
//
// dim is the number of floats per point; semantic says whether coordinate conversions apply. The
// name is a free-form label (e.g. "roughness", or a PLY property name) used to look it up.
//
// Preconditions for saving (checked; saveSpz fails and logs the problem otherwise):
//   - name is 1-65535 bytes of printable ASCII (0x20-0x7E), and unique among the cloud's custom
//     attributes.
//   - dim >= 1, and a multiple of 3 for semantic Spatial.
//   - bits is 32 (float32, lossless) or 1-16 (each component quantized to 2^bits levels between
//     its min and max over all points; values must then be finite).
//   - values holds numPoints * dim floats, point-major.
//   - The file version is 4 or later.
struct SpzExtensionCustomAttributeNiantic : public SpzExtensionBase,
                                            public SpzFramedExtension,
                                            public SpzKeyedExtension,
                                            public SpzCoordinateDependentExtension {
  std::string name;
  uint16_t dim = 1;
  SpzAttributeSemantic semantic = SpzAttributeSemantic::Custom;
  uint8_t bits = 32;
  std::vector<float> values;

  SpzExtensionCustomAttributeNiantic();

  // A custom attribute's values only fit in extension frames, which saveSpz writes through
  // planFrameWrites(). write(os) therefore writes nothing (and logs an error), and payloadBytes()
  // is 0.
  uint32_t payloadBytes() const override;
  void write(std::ostream& os) const override;
  SpzExtensionBase* copyAsRawData() const override;

  // PLY has no representation for custom attributes yet: they aren't read from PLY, and saving to
  // PLY drops them with a warning.
  std::optional<std::shared_ptr<SpzExtensionBase>> tryReadFromPly(
      std::istream& in, const std::vector<PlyExtraElement>& elements) const override;
  void writePlyHeader(std::ostream& out) const override;
  void writePlyData(std::ostream& out) const override;

  std::optional<SpzFrameWritePlan> planFrameWrites(int32_t numPoints) const override;
  std::string uniquenessKey() const override;
  SpzExtensionBasePtr withConvertedCoordinates(const CoordinateConverter& c) const override;

  // This extension's type ID, SPZ_NIANTIC_custom_attribute.
  static SpzExtensionType type();
};

// Returns the custom attribute called name (compared byte for byte), or nullptr if extensions has
// none. Coordinate conversions replace SPATIAL attributes with converted copies, so look the
// attribute up again after converting the cloud.
std::shared_ptr<SpzExtensionCustomAttributeNiantic> findCustomAttribute(
  const std::vector<SpzExtensionBasePtr>& extensions, const std::string& name);

}  // namespace spz

#endif  // SPZ_CUSTOM_ATTRIBUTE_NIANTIC_H_
