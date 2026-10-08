/*
MIT License

Copyright (c) 2026 Niantic Labs

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

#ifndef SPZ_GEOREFERENCE_NIANTIC_H_
#define SPZ_GEOREFERENCE_NIANTIC_H_

#include <array>
#include <limits>
#include <string>

#include "splat-extensions.h"

namespace spz {

// Records a similarity transform that places the asset in a geocentric CRS:
//
//   p_crs = scale * (q * p_local * q^-1) + origin
//
// p_local is a position in RUB, whatever frame the data is stored in. q is `rotation`, a unit
// Hamilton quaternion (x, y, z, w) applied as an active rotation; `scale` is uniform and positive;
// `origin` is in meters. `crs` names the target: "AUTHORITY:CODE" when the target has a code,
// otherwise an inline PROJJSON GeodeticCRS. See extensions/README.md for the full rules.
//
// This extension is a descriptor: the library never applies the transform to Gaussian data.
struct SpzExtensionGeoreferenceNiantic : public SpzExtensionBase {
  // Serialization of `crs`, stored in flags bits 1-2. Unassigned values reject the payload.
  enum class CrsEncoding : uint8_t { AuthorityCode = 0, Projjson = 1 };

  static constexpr uint8_t kExtVersion = 1;
  static constexpr uint8_t kFlagHasEpoch = 1 << 0;
  static constexpr uint8_t kCrsEncodingMask = 0x06;  // bits 1-2
  static constexpr uint8_t kCrsEncodingShift = 1;
  static constexpr uint16_t kMaxCrsBytes = 0xFFFF;  // u16 crs_length

  std::string crs;  // Target CRS, per crsEncoding; required
  CrsEncoding crsEncoding = CrsEncoding::AuthorityCode;
  std::array<double, 3> origin = {0.0, 0.0, 0.0};          // Meters
  std::array<double, 4> rotation = {0.0, 0.0, 0.0, 1.0};   // Unit quaternion x, y, z, w
  double scale = 1.0;                                      // Finite and > 0
  // Coordinate epoch of p_crs as a decimal year; NaN when absent.
  double epoch = std::numeric_limits<double>::quiet_NaN();

  SpzExtensionGeoreferenceNiantic();
  uint32_t payloadBytes() const override;
  void write(std::ostream& os) const override;
  SpzExtensionBase* copyAsRawData() const override;
  std::optional<std::shared_ptr<SpzExtensionBase>> tryReadFromPly(
      std::istream& in, const std::unordered_set<std::string>& elementNames) const override;
  void writePlyHeader(std::ostream& out) const override;
  void writePlyData(std::ostream& out) const override;
  static std::optional<SpzExtensionBasePtr> read(std::istream& is);
  static SpzExtensionType type();
};

}  // namespace spz

#endif  // SPZ_GEOREFERENCE_NIANTIC_H_
