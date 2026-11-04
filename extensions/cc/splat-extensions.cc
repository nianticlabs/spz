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
#include "georeference-niantic.h"
#include "custom-attribute-niantic.h"
#include "extension-frames.h"
#include "environment-map-adobe.h"
#include "load-spz.h"

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_set>

namespace spz {
std::streamoff detail::remainingBytes(std::istream& is) {
  const std::streampos cur = is.tellg();
  if (cur == std::streampos(-1)) return -1;
  is.seekg(0, std::ios::end);
  const std::streampos end = is.tellg();
  is.seekg(cur);
  if (end == std::streampos(-1)) return -1;
  return end - cur;
}

namespace {
using detail::readExact;
using detail::remainingBytes;

// Extension stream format: [u32 type][u32 byteLength][payload...] per record, until EOF.
// Unknown types are skipped by advancing byteLength bytes so multiple vendors can coexist.
bool tryParseExtension(std::istream& is, std::vector<SpzParsedExtension>& out) {
  uint32_t type_u32{};
  if (!readExact(is, type_u32))
    return false;  // EOF or short read -> stop
  uint32_t byteLength{};
  if (!readExact(is, byteLength))
    return false;

  // Bound byteLength against bytes actually remaining in the stream before any allocation,
  // so a malformed file cannot drive a multi-GB allocation that the subsequent read would reject.
  const std::streamoff remaining = remainingBytes(is);
  if (remaining < 0 || static_cast<uint64_t>(byteLength) > static_cast<uint64_t>(remaining)) {
    SpzLog("[SPZ ERROR] tryParseExtension: byteLength %u exceeds remaining stream bytes",
           static_cast<unsigned>(byteLength));
    return false;
  }

  SpzExtensionType type = static_cast<SpzExtensionType>(type_u32);

  switch (type) {
    case SpzExtensionType::SPZ_ADOBE_safe_orbit_camera: {
      std::vector<char> payload(byteLength);
      if (!is.read(payload.data(), static_cast<std::streamsize>(byteLength)))
        return false;
      std::istringstream iss(std::string(payload.data(), payload.size()));
      auto rec = SpzExtensionSafeOrbitCameraAdobe::read(iss);
      if (rec)
        out.push_back({std::move(*rec), std::nullopt});
      return true;  // continue to next record
    }
    case SpzExtensionType::SPZ_ADOBE_coordinate_system: {
      std::vector<char> payload(byteLength);
      if (!is.read(payload.data(), static_cast<std::streamsize>(byteLength)))
        return false;
      std::istringstream iss(std::string(payload.data(), payload.size()));
      auto rec = SpzExtensionCoordinateSystemAdobe::read(iss);
      if (rec)
        out.push_back({std::move(*rec), std::nullopt});
      return true;
    }
    case SpzExtensionType::SPZ_NIANTIC_custom_attribute: {
      std::vector<char> payload(byteLength);
      if (!is.read(payload.data(), static_cast<std::streamsize>(byteLength)))
        return false;
      std::istringstream iss(std::string(payload.data(), payload.size()));
      auto readPlan = parseCustomAttributeNiantic(iss);
      if (readPlan)
        out.push_back({nullptr, std::move(*readPlan)});
      return true;
    }
    case SpzExtensionType::SPZ_NIANTIC_georeference: {
      std::vector<char> payload(byteLength);
      if (!is.read(payload.data(), static_cast<std::streamsize>(byteLength)))
        return false;
      std::istringstream iss(std::string(payload.data(), payload.size()));
      auto rec = SpzExtensionGeoreferenceNiantic::read(iss);
      if (rec)
        out.push_back({std::move(*rec), std::nullopt});
      return true;
    }
    case SpzExtensionType::SPZ_ADOBE_environment_map: {
      std::string payload(byteLength, '\0');
      if (!is.read(payload.data(), static_cast<std::streamsize>(byteLength)))
        return false;
      std::istringstream iss(payload);
      auto rec = SpzExtensionEnvironmentMapAdobe::read(iss);
      if (rec)
        out.push_back({std::move(*rec), std::nullopt});
      return true;
    }
    default:
      // Unknown type: skip payload and continue
      SpzLog("[SPZ WARNING] Unknown extension type 0x%08x (%u bytes) was skipped — loaded data may be incorrect",
             static_cast<unsigned>(type_u32), static_cast<unsigned>(byteLength));
      if (byteLength > 0)
        is.ignore(static_cast<std::streamsize>(byteLength));
      return true;
  }
}

using PlyExemplarEntry = std::pair<const std::unordered_set<std::string>*, SpzExtensionBasePtr>;

static const std::vector<PlyExemplarEntry>& getPlyExtensionRegistry() {
  static const std::vector<PlyExemplarEntry> registry = {
    {&SpzExtensionSafeOrbitCameraAdobe::kRequiredPlyElementNames, std::make_shared<SpzExtensionSafeOrbitCameraAdobe>()},
    {&SpzExtensionEnvironmentMapAdobe::kRequiredPlyElementNames, std::make_shared<SpzExtensionEnvironmentMapAdobe>()},
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

CoordinateSystem getPackedCoordinateSystem(
    const std::vector<SpzExtensionBasePtr>& extensions) {
  if (auto coordExt = findExtensionByType<SpzExtensionCoordinateSystemAdobe>(extensions))
    return coordExt->resolve();
  // Unrecognized extension types are already warned about in tryParseExtension.
  return CoordinateSystem::RUB;  // default: no known extension overrides the packed coordinate system
}

bool allExtensionsValid(const std::vector<SpzExtensionBasePtr>& extensions) {
  for (const auto& ext : extensions) {
    if (!ext->isValid()) return false;
  }
  return true;
}

bool isKnownPlyExtensionElement(const std::string& elementName) {
  return getKnownPlyExtensionElementNames().count(elementName) != 0;
}

void readExtensionRecords(std::istream& is, std::vector<SpzParsedExtension>& out) {
  while (tryParseExtension(is, out));

  // Stopped because there are no more extension records (EOF). Clear stream state
  // so the caller does not treat expected end-of-extensions as a read error.
  if (is.eof()) {
    is.clear();
  }
}

void readAllExtensions(std::istream& is, std::vector<SpzExtensionBasePtr>& out) {
  std::vector<SpzParsedExtension> parsed;
  readExtensionRecords(is, parsed);
  for (auto& entry : parsed) {
    if (entry.record) {
      out.push_back(std::move(entry.record));
    } else if (entry.readPlan) {
      SpzLog("[SPZ ERROR] Dropped extension 0x%08x: its data is only stored in SPZ version 4 files",
             static_cast<unsigned>(entry.readPlan->type));
    }
  }
}

void writeAllExtensions(const std::vector<SpzExtensionBasePtr>& list, std::ostream& os) {
  for (const auto& rec : list) {
    if (!rec) continue;
    if (dynamic_cast<const SpzFramedExtension*>(rec.get())) {
      SpzLog("[SPZ ERROR] Skipped extension 0x%08x: its data can only be saved by saveSpz in SPZ "
             "version 4 files", static_cast<unsigned>(rec->extensionType));
      continue;
    }
    rec->write(os);
  }
}

void readExtensionsFromPly(std::istream& in, const std::vector<spz::PlyExtraElement>& extraElements, std::vector<SpzExtensionBasePtr>& extensions) {
  const auto& registry = getPlyExtensionRegistry();
  for (auto it = extraElements.begin(); it != extraElements.end();) {
    const auto entry = std::find_if(registry.begin(), registry.end(),
                                    [&](const auto& e) { return e.first->count(it->name) != 0; });
    // A known element starts the run of consecutive elements its extension declares.
    const auto end = entry == registry.end()
        ? std::next(it)
        : std::find_if(it, extraElements.end(), [&](const auto& e) { return entry->first->count(e.name) == 0; });
    uint64_t bytes = 0;
    for (auto e = it; e != end; ++e)
      bytes += static_cast<uint64_t>(e->count) * e->bytesPerElement;
    const std::streampos start = in.tellg();
    if (entry != registry.end()) {
      auto rec = entry->second->tryReadFromPly(in, std::vector<PlyExtraElement>(it, end));
      if (rec)
        extensions.push_back(std::move(*rec));
    } else {
      SpzLog("[SPZ] Skipped %llu bytes for element '%s'", static_cast<unsigned long long>(bytes), it->name.c_str());
    }
    // Continue after the run whether or not the extension consumed all of it.
    in.clear();
    in.seekg(start + static_cast<std::streamoff>(bytes));
    it = end;
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
