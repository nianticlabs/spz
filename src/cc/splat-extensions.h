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
  float shMin = -1.0f;     // Minimum SH coefficient value used for quantization
  float shMax = 1.0f;      // Maximum SH coefficient value used for quantization

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