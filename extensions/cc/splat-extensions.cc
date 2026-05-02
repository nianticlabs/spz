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
#include "safe-orbit-camera-adobe.h"
#include "coordinate-system-adobe.h"
#include "load-spz.h"

#include <sstream>
#include <string>
#include <unordered_set>

namespace spz {
namespace {
template <class T>
bool readExact(std::istream& is, T& out) {
  return static_cast<bool>(is.read(reinterpret_cast<char*>(&out), sizeof(T)));
}

// Extension stream format: [u32 type][u32 byteLength][payload...] per record, until EOF.
// Unknown types are skipped by advancing byteLength bytes so multiple vendors can coexist.
bool tryParseExtension(std::istream& is, std::vector<SpzExtensionBasePtr>& out) {
  uint32_t type_u32{};
  if (!readExact(is, type_u32))
    return false;  // EOF or short read -> stop
  uint32_t byteLength{};
  if (!readExact(is, byteLength))
    return false;

  SpzExtensionType type = static_cast<SpzExtensionType>(type_u32);

  switch (type) {
    case SpzExtensionType::SPZ_ADOBE_safe_orbit_camera: {
      std::vector<char> payload(byteLength);
      if (!is.read(payload.data(), static_cast<std::streamsize>(byteLength)))
        return false;
      std::istringstream iss(std::string(payload.data(), payload.size()));
      auto rec = SpzExtensionSafeOrbitCameraAdobe::read(iss);
      if (rec)
        out.push_back(std::move(*rec));
      return true;  // continue to next record
    }
    case SpzExtensionType::SPZ_ADOBE_coordinate_system: {
      std::vector<char> payload(byteLength);
      if (!is.read(payload.data(), static_cast<std::streamsize>(byteLength)))
        return false;
      std::istringstream iss(std::string(payload.data(), payload.size()));
      auto rec = SpzExtensionCoordinateSystemAdobe::read(iss);
      if (rec)
        out.push_back(std::move(*rec));
      return true;
    }
    default:
      // Unknown type: skip payload and continue
      SpzLog("[SPZ WARNING] Skipping unknown extension type: 0x%08x (%u bytes)", static_cast<unsigned>(type_u32), static_cast<unsigned>(byteLength));
      if (byteLength > 0)
        is.ignore(static_cast<std::streamsize>(byteLength));
      return true;
  }
}

using PlyExemplarEntry = std::pair<const std::unordered_set<std::string>*, SpzExtensionBasePtr>;

static const std::vector<PlyExemplarEntry>& getPlyExtensionRegistry() {
  static const std::vector<PlyExemplarEntry> registry = {
    {&SpzExtensionSafeOrbitCameraAdobe::kRequiredPlyElementNames, std::make_shared<SpzExtensionSafeOrbitCameraAdobe>()},
  };
  return registry;
}

const std::unordered_set<std::string>& getKnownPlyExtensionElementNames() {
  static const std::unordered_set<std::string> known = []() {
    std::unordered_set<std::string> out;
    for (const auto& entry : getPlyExtensionRegistry())
      for (const auto& name : *entry.first)
        out.insert(name);
    return out;
  }();
  return known;
}
}  // namespace

bool isKnownPlyExtensionElement(const std::string& elementName) {
  return getKnownPlyExtensionElementNames().count(elementName) != 0;
}

void readAllExtensions(std::istream& is, std::vector<SpzExtensionBasePtr>& out) {
  while (tryParseExtension(is, out));

  // Stopped because there are no more extension records (EOF). Clear stream state
  // so the caller does not treat expected end-of-extensions as a read error.
  if (is.eof()) {
    is.clear();
  }
}

void writeAllExtensions(const std::vector<SpzExtensionBasePtr>& list, std::ostream& os) {
  for (const auto& rec : list)
    rec->write(os);
}

void readExtensionsFromPly(std::istream& in, const std::vector<spz::PlyExtraElement>& extraElements, std::vector<SpzExtensionBasePtr>& extensions) {
  std::unordered_set<std::string> elementNames;
  for (const auto& elem : extraElements)
    elementNames.insert(elem.name);

  for (const auto& entry : getPlyExtensionRegistry()) {
    auto rec = entry.second->tryReadFromPly(in, elementNames);
    if (rec)
      extensions.push_back(std::move(*rec));
  }
}

void writeExtensionsToPlyHeader(const std::vector<SpzExtensionBasePtr>& extensions, std::ostream& out) {
  for (const auto& ext : extensions)
    ext->writePlyHeader(out);
}

void writeExtensionsToPlyData(const std::vector<SpzExtensionBasePtr>& extensions, std::ostream& out) {
  for (const auto& ext : extensions)
    ext->writePlyData(out);
}
}  // namespace spz
