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

#ifndef SPZ_UNITS_ADOBE_H_
#define SPZ_UNITS_ADOBE_H_

#include "splat-extensions.h"

namespace spz {

// A length unit recorded by SPZ_ADOBE_units.
enum class LengthUnit : uint32_t {
  UNSPECIFIED = 0,  // Not valid in a record
  MILLIMETER = 1,
  CENTIMETER = 2,
  METER = 3,
  KILOMETER = 4,
};

// Records the length unit the positions and scales are written in. A label only: the library never
// converts the data, and other extensions don't refer to it. Clients that care scale the loaded
// data themselves (by 0.001 for millimeters to meters). Without this record the unit is unspecified;
// by SPZ convention, meters.
struct SpzExtensionUnitsAdobe : public SpzExtensionBase {
  LengthUnit lengthUnit = LengthUnit::UNSPECIFIED;

  SpzExtensionUnitsAdobe();
  uint32_t payloadBytes() const override;
  void write(std::ostream& os) const override;
  // False (after logging why) unless lengthUnit is MILLIMETER, CENTIMETER, METER or KILOMETER.
  bool isValid() const override;
  SpzExtensionBase* copyAsRawData() const override;
  std::optional<std::shared_ptr<SpzExtensionBase>> tryReadFromPly(
      std::istream& in, const std::vector<PlyExtraElement>& elements) const override;
  void writePlyHeader(std::ostream& out) const override;
  void writePlyData(std::ostream& out) const override;
  static std::optional<SpzExtensionBasePtr> read(std::istream& is);
  static SpzExtensionType type();
};

}  // namespace spz

#endif  // SPZ_UNITS_ADOBE_H_
