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
#include "splat-types.h"

namespace spz {

// Records the length unit in which the Gaussian data is stored in the file, e.g. millimeters for a
// finer position grid and smaller representable scales.
//
// Like SPZ_ADOBE_coordinate_system, this is a descriptor, not an instruction. packGaussians converts
// PackOptions::fromUnit to lengthUnit (instead of the default meters), and unpackGaussians converts
// lengthUnit to UnpackOptions::toUnit (meters by default). Consumers using those functions must not
// scale the data again based on this extension.
//
// UNSPECIFIED is equivalent to the extension being absent: both mean meters, but a present
// UNSPECIFIED value is likely a writer bug and will produce a warning.
struct SpzExtensionUnitsAdobe : public SpzExtensionBase {
  LengthUnit lengthUnit = LengthUnit::UNSPECIFIED;

  SpzExtensionUnitsAdobe();
  uint32_t payloadBytes() const override;
  void write(std::ostream& os) const override;
  // False (after logging why) if lengthUnit is not a LengthUnit value.
  bool isValid() const override;
  SpzExtensionBase* copyAsRawData() const override;
  std::optional<std::shared_ptr<SpzExtensionBase>> tryReadFromPly(
      std::istream& in, const std::vector<PlyExtraElement>& elements) const override;
  void writePlyHeader(std::ostream& out) const override;
  void writePlyData(std::ostream& out) const override;
  static std::optional<SpzExtensionBasePtr> read(std::istream& is);
  static SpzExtensionType type();

  // Returns the storage length unit recorded by this extension.
  // Returns METER and logs a warning if lengthUnit is UNSPECIFIED (likely a writer bug).
  LengthUnit resolve() const;
};

}  // namespace spz

#endif  // SPZ_UNITS_ADOBE_H_
