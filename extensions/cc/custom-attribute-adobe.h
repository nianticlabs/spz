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

#ifndef SPZ_CUSTOM_ATTRIBUTE_ADOBE_H_
#define SPZ_CUSTOM_ATTRIBUTE_ADOBE_H_

#include "splat-extensions.h"

namespace spz {

// How a custom attribute's values follow coordinate-system conversions.
enum class SpzCustomAttributeKind : uint8_t {
  Scalar = 0,   // never converted
  Vector3 = 1,  // each (x, y, z) triple is a direction, converted like positions
};

// One named per-point attribute outside the core set (normals, material parameters, feature
// vectors, ...). Files carry one record per attribute, any number of them.
//
// `values` holds numPoints * dim floats, point-major. `bits` sets the stored precision: 32 keeps
// float32 exactly; 1-8 quantizes each component to 2^bits levels between its min and max over all
// points. The values are zstd-compressed in the record.
struct SpzExtensionCustomAttributeAdobe : public SpzExtensionBase {
  std::string name;
  uint16_t dim = 1;
  SpzCustomAttributeKind kind = SpzCustomAttributeKind::Scalar;
  uint8_t bits = 32;
  std::vector<float> values;

  SpzExtensionCustomAttributeAdobe();
  uint32_t payloadBytes() const override;
  void write(std::ostream& os) const override;
  SpzExtensionBase* copyAsRawData() const override;
  std::optional<std::shared_ptr<SpzExtensionBase>> tryReadFromPly(
      std::istream& in, const std::vector<PlyExtraElement>& elements) const override;
  void writePlyHeader(std::ostream& out) const override;
  void writePlyData(std::ostream& out) const override;
  // numPoints comes from the file header; records that don't cover exactly that many points are rejected.
  static std::optional<SpzExtensionBasePtr> read(std::istream& is, int32_t numPoints);
  static SpzExtensionType type();

  // Returns false and logs the problem if this attribute can't be saved with a cloud of numPoints.
  bool validate(int32_t numPoints) const;

 private:
  std::string encode() const;
};

}  // namespace spz

#endif  // SPZ_CUSTOM_ATTRIBUTE_ADOBE_H_
