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

#include "units-adobe.h"
#include "load-spz.h"

namespace spz {
namespace {
bool isKnownLengthUnit(uint32_t value) {
  return value >= static_cast<uint32_t>(LengthUnit::MILLIMETER) &&
         value <= static_cast<uint32_t>(LengthUnit::KILOMETER);
}
}  // namespace

SpzExtensionUnitsAdobe::SpzExtensionUnitsAdobe()
    : SpzExtensionBase(SpzExtensionType::SPZ_ADOBE_units) {}

uint32_t SpzExtensionUnitsAdobe::payloadBytes() const {
  return sizeof(uint32_t);  // LengthUnit value
}

void SpzExtensionUnitsAdobe::write(std::ostream& os) const {
  const uint32_t t = static_cast<uint32_t>(extensionType);
  const uint32_t len = payloadBytes();
  const uint32_t unitValue = static_cast<uint32_t>(lengthUnit);
  SpzLog("[SPZ] Writing extension: UnitsAdobe (unit=%u)", static_cast<unsigned>(unitValue));
  os.write(reinterpret_cast<const char*>(&t), sizeof(t));
  os.write(reinterpret_cast<const char*>(&len), sizeof(len));
  os.write(reinterpret_cast<const char*>(&unitValue), sizeof(unitValue));
}

bool SpzExtensionUnitsAdobe::isValid() const {
  if (!isKnownLengthUnit(static_cast<uint32_t>(lengthUnit))) {
    SpzLog("[SPZ ERROR] UnitsAdobe: length unit must be MILLIMETER, CENTIMETER, METER or KILOMETER, got %u",
           static_cast<unsigned>(lengthUnit));
    return false;
  }
  return true;
}

std::optional<SpzExtensionBasePtr> SpzExtensionUnitsAdobe::read(std::istream& is) {
  uint32_t unitValue{};
  if (!detail::readExact(is, unitValue)) {
    SpzLog("[SPZ WARNING] UnitsAdobe: payload too short, record dropped");
    return std::nullopt;
  }
  if (!isKnownLengthUnit(unitValue)) {
    SpzLog("[SPZ WARNING] UnitsAdobe: unknown length unit %u, record dropped", static_cast<unsigned>(unitValue));
    return std::nullopt;
  }
  auto rec = std::make_shared<SpzExtensionUnitsAdobe>();
  rec->lengthUnit = static_cast<LengthUnit>(unitValue);
  SpzLog("[SPZ] Found extension: UnitsAdobe (unit=%u)", static_cast<unsigned>(unitValue));
  return std::optional{std::move(rec)};
}

SpzExtensionType SpzExtensionUnitsAdobe::type() {
  return SpzExtensionType::SPZ_ADOBE_units;
}

SpzExtensionBase* SpzExtensionUnitsAdobe::copyAsRawData() const {
  return new SpzExtensionUnitsAdobe(*this);
}

// No PLY representation for this extension.
std::optional<SpzExtensionBasePtr> SpzExtensionUnitsAdobe::tryReadFromPly(
    std::istream&, const std::vector<PlyExtraElement>&) const {
  return std::nullopt;
}

void SpzExtensionUnitsAdobe::writePlyHeader(std::ostream&) const {}

void SpzExtensionUnitsAdobe::writePlyData(std::ostream&) const {}

}  // namespace spz
