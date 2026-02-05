/*
MIT License

Copyright (c) 2025 Niantic Labs
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

#ifndef SPZ_SPLAT_EXTENSIONS_H_
#define SPZ_SPLAT_EXTENSIONS_H_

#include <cstdint>
#include <iostream>
#include <istream>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <vector>

namespace spz {
enum class SpzExtensionType : uint32_t {
  SPZ_ADOBE_sh_quantization = 0xADBE0001u,
  SPZ_ADOBE_safe_orbit_camera = 0xADBE0002u,
};

struct SpzExtensionBase {
  SpzExtensionType extensionType;
  explicit SpzExtensionBase(SpzExtensionType t) : extensionType(t) {
  }
  virtual ~SpzExtensionBase() = default;
  virtual void write(std::ostream& os) const = 0;
  virtual SpzExtensionBase* copyAsRawData() const = 0;
};

using SpzExtensionBasePtr = std::shared_ptr<SpzExtensionBase>;

// SH quantization parameters
struct SpzExtensionSHQuantizationAdobe : public SpzExtensionBase {
  uint8_t sh1Bits = 5;     // Bits for SH degree 1 coefficients
  uint8_t shRestBits = 4;  // Bits for SH degree 2+ coefficients

  SpzExtensionSHQuantizationAdobe();
  void write(std::ostream& os) const override;
  SpzExtensionBase* copyAsRawData() const override;
  static std::optional<SpzExtensionBasePtr> read(std::istream& is);
  static SpzExtensionType type();
};

struct SpzExtensionSafeOrbitCameraAdobe : public SpzExtensionBase {
  float safeOrbitElevationMin = 0.0f;  // Minimum elevation for safe orbit (radians)
  float safeOrbitElevationMax = 0.0f;  // Maximum elevation for safe orbit (radians)
  float safeOrbitRadiusMin = 0.0f;     // Minimum radius for safe orbit

  SpzExtensionSafeOrbitCameraAdobe();
  void write(std::ostream& os) const override;
  SpzExtensionBase* copyAsRawData() const override;
  static std::optional<SpzExtensionBasePtr> read(std::istream& is);
  static SpzExtensionType type();
};

void readAllExtensions(std::istream& is, std::vector<SpzExtensionBasePtr>& out);

void writeAllExtensions(const std::vector<SpzExtensionBasePtr>& list, std::ostream& os);

template <typename T>
std::shared_ptr<T> findExtensionByType(const std::vector<SpzExtensionBasePtr>& list) {
  for (const auto& rec : list) {
    if (rec->extensionType == T::type())
      return std::dynamic_pointer_cast<T>(rec);
  }
  return nullptr;
}
}  // namespace spz

#endif  // SPZ_SPLAT_EXTENSIONS_H_