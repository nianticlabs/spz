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

#include "splat-extensions.h"

#include "load-spz.h"

namespace spz {
namespace {
template <class T>
bool readExact(std::istream& is, T& out) {
  return static_cast<bool>(is.read(reinterpret_cast<char*>(&out), sizeof(T)));
}

struct StreamMark {
  std::istream& is;
  std::streampos pos;
  bool committed = false;
  explicit StreamMark(std::istream& s) : is(s), pos(s.tellg()) {
  }
  void commit() {
    committed = true;
  }
  ~StreamMark() {
    if (!committed) {
      is.clear();
      is.seekg(pos);
    }
  }
};

std::optional<SpzExtensionBasePtr> tryParseExtension(std::istream& is) {
  StreamMark mark(is);

  uint32_t type_u32{};
  if (!readExact(is, type_u32))
    return std::nullopt;  // EOF/short
  SpzExtensionType type = static_cast<SpzExtensionType>(type_u32);

  switch (type) {
    case SpzExtensionType::SPZ_ADOBE_sh_quantization:
      mark.commit();
      return SpzExtensionSHQuantizationAdobe::read(is);
    case SpzExtensionType::SPZ_ADOBE_safe_orbit_camera:
      mark.commit();
      return SpzExtensionSafeOrbitCameraAdobe::read(is);
    default:
      // Unknown tag -> not one of our records.
      return std::nullopt;
  }
}
}  // namespace

SpzExtensionSHQuantizationAdobe::SpzExtensionSHQuantizationAdobe() : SpzExtensionBase(SpzExtensionType::SPZ_ADOBE_sh_quantization) {
}

void SpzExtensionSHQuantizationAdobe::write(std::ostream& os) const {
  const uint32_t t = static_cast<uint32_t>(extensionType);
  SpzLog("[SPZ] Writing extension: SHQuantizationAdobe");
  os.write(reinterpret_cast<const char*>(&t), sizeof(t));
  os.write(reinterpret_cast<const char*>(&sh1Bits), sizeof(sh1Bits));
  os.write(reinterpret_cast<const char*>(&shRestBits), sizeof(shRestBits));
}

SpzExtensionBase* SpzExtensionSHQuantizationAdobe::copyAsRawData() const {
  return new SpzExtensionSHQuantizationAdobe(*this);
}

std::optional<SpzExtensionBasePtr> SpzExtensionSHQuantizationAdobe::read(std::istream& is) {
  StreamMark mark(is);
  SpzLog("[SPZ] Found extension: SpzExtensionSHQuantizationAdobe");
  auto rec = std::make_shared<SpzExtensionSHQuantizationAdobe>();
  if (!readExact(is, rec->sh1Bits) || !readExact(is, rec->shRestBits)) {
    SpzLog("[SPZ ERROR] Failed to read all fields for SpzExtensionSHQuantizationAdobe");
    return std::nullopt;
  }

  mark.commit();
  return std::optional{std::move(rec)};
}

SpzExtensionType SpzExtensionSHQuantizationAdobe::type() {
  return SpzExtensionType::SPZ_ADOBE_sh_quantization;
}

SpzExtensionSafeOrbitCameraAdobe::SpzExtensionSafeOrbitCameraAdobe() : SpzExtensionBase(SpzExtensionType::SPZ_ADOBE_safe_orbit_camera) {
}

void SpzExtensionSafeOrbitCameraAdobe::write(std::ostream& os) const {
  const uint32_t t = static_cast<uint32_t>(extensionType);
  SpzLog("[SPZ] Writing extension: SafeOrbitCameraAdobe");
  os.write(reinterpret_cast<const char*>(&t), sizeof(t));
  os.write(reinterpret_cast<const char*>(&safeOrbitElevationMin), sizeof(safeOrbitElevationMin));
  os.write(reinterpret_cast<const char*>(&safeOrbitElevationMax), sizeof(safeOrbitElevationMax));
  os.write(reinterpret_cast<const char*>(&safeOrbitRadiusMin), sizeof(safeOrbitRadiusMin));
}

std::optional<SpzExtensionBasePtr> SpzExtensionSafeOrbitCameraAdobe::read(std::istream& is) {
  StreamMark mark(is);
  SpzLog("[SPZ] Found extension: SafeOrbitCameraAdobe");
  auto rec = std::make_shared<SpzExtensionSafeOrbitCameraAdobe>();
  if (!readExact(is, rec->safeOrbitElevationMin) || !readExact(is, rec->safeOrbitElevationMax) || !readExact(is, rec->safeOrbitRadiusMin)) {
    SpzLog("[SPZ ERROR] Failed to read all fields for SafeOrbitCameraAdobe");
    return std::nullopt;
  }

  mark.commit();
  return std::optional{std::move(rec)};
}

SpzExtensionType SpzExtensionSafeOrbitCameraAdobe::type() {
  return SpzExtensionType::SPZ_ADOBE_safe_orbit_camera;
}

SpzExtensionBase* SpzExtensionSafeOrbitCameraAdobe::copyAsRawData() const {
  return new SpzExtensionSafeOrbitCameraAdobe(*this);
}

void readAllExtensions(std::istream& is, std::vector<SpzExtensionBasePtr>& out) {
  while (true) {
    auto rec = tryParseExtension(is);
    if (!rec)
      break;
    out.push_back(std::move(*rec));
  }
}

void writeAllExtensions(const std::vector<SpzExtensionBasePtr>& list, std::ostream& os) {
  for (const auto& rec : list)
    rec->write(os);
}
}  // namespace spz
