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

#include "environment-map-adobe.h"
#include "extension-frames.h"
#include "load-spz.h"
#include "splat-types.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>

namespace spz {
namespace {
constexpr float kMaxHalf = 65504.0f;
// width, height, origin, up, e1, e2, metadataLength
constexpr size_t kFixedBytes = 2 * sizeof(uint32_t) + 12 * sizeof(float) + sizeof(uint32_t);
constexpr uint32_t kChannels = 3;  // RGB
// PLY elements, written in this order after the vertex element.
constexpr const char* kPlyEnvironment = "environment";
constexpr const char* kPlyMap = "environment_map";
constexpr const char* kPlyMetadata = "environment_metadata";

template <class T>
bool readExact(std::istream& is, T& out) {
  return static_cast<bool>(is.read(reinterpret_cast<char*>(&out), sizeof(T)));
}

uint64_t pixelFloats(uint32_t width, uint32_t height) {
  return static_cast<uint64_t>(width) * height * kChannels;
}

// Byte size of the float16 pixels, or std::nullopt if it can't be addressed on this platform.
std::optional<size_t> pixelBytes(uint32_t width, uint32_t height) {
  const uint64_t bytes = pixelFloats(width, height) * sizeof(Half);
  if (bytes > std::numeric_limits<size_t>::max()) return std::nullopt;
  return static_cast<size_t>(bytes);
}

// Everything in the record but the frame table: width, height, the frame vectors and the metadata.
std::string encodeHeader(const SpzExtensionEnvironmentMapAdobe& env) {
  std::ostringstream os;
  writeValue(os, env.width);
  writeValue(os, env.height);
  for (const auto* v : {&env.origin, &env.up, &env.e1, &env.e2})
    os.write(reinterpret_cast<const char*>(v->data()), v->size() * sizeof(float));
  writeValue(os, static_cast<uint32_t>(env.metadata.size()));
  os.write(env.metadata.data(), static_cast<std::streamsize>(env.metadata.size()));
  return os.str();
}

bool isFinite(const std::array<float, 3>& v) {
  return std::all_of(v.begin(), v.end(), [](float x) { return std::isfinite(x); });
}
}  // namespace

const std::unordered_set<std::string> SpzExtensionEnvironmentMapAdobe::kRequiredPlyElementNames = {
  kPlyEnvironment,
  kPlyMap,
  kPlyMetadata,
};

SpzExtensionEnvironmentMapAdobe::SpzExtensionEnvironmentMapAdobe()
    : SpzExtensionBase(SpzExtensionType::SPZ_ADOBE_environment_map) {}

bool SpzExtensionEnvironmentMapAdobe::isValid() const {
  const char* problem = nullptr;
  if (width == 0 || height == 0)
    problem = "width and height must be at least 1";
  else if (pixels.size() != pixelFloats(width, height))
    problem = "pixels must hold width * height * 3 floats";
  else if (!isFinite(origin) || !isFinite(up) || !isFinite(e1) || !isFinite(e2))
    problem = "origin, up, e1 and e2 must be finite";
  else if (!std::all_of(pixels.begin(), pixels.end(),
                        [](float v) { return std::isfinite(v) && std::abs(v) <= kMaxHalf; }))
    problem = "pixels must be finite and within the float16 range";
  else if (kFixedBytes + metadata.size() > std::numeric_limits<uint32_t>::max() / 2)
    problem = "metadata is too large for an extension record";
  if (problem) {
    SpzLog("[SPZ ERROR] Environment map: %s", problem);
    return false;
  }
  return true;
}

std::string SpzExtensionEnvironmentMapAdobe::uniquenessKey() const { return "environment_map"; }

SpzExtensionBasePtr SpzExtensionEnvironmentMapAdobe::withConvertedCoordinates(const CoordinateConverter& c) const {
  if (isIdentityPositionTransform(c)) return nullptr;
  auto converted = std::make_shared<SpzExtensionEnvironmentMapAdobe>(*this);
  for (auto* v : {&converted->origin, &converted->up, &converted->e1, &converted->e2})
    convertDirection(c, v->data());
  return converted;
}

uint32_t SpzExtensionEnvironmentMapAdobe::payloadBytes() const { return 0; }

void SpzExtensionEnvironmentMapAdobe::write(std::ostream& /*os*/) const {
  SpzLog("[SPZ ERROR] Environment map not written: its pixels can only be saved by saveSpz in SPZ version 4 files");
}

std::optional<SpzFrameWritePlan> SpzExtensionEnvironmentMapAdobe::planFrameWrites(int32_t /*numPoints*/) const {
  const std::optional<size_t> bytes = pixelBytes(width, height);
  if (!isValid() || !bytes) return std::nullopt;
  std::vector<uint8_t> data(*bytes);
  for (size_t i = 0; i < pixels.size(); ++i) {
    const Half half = floatToHalf(pixels[i]);
    std::memcpy(&data[i * sizeof(Half)], &half, sizeof(Half));
  }
  SpzLog("[SPZ] Writing extension: EnvironmentMapAdobe (%ux%u)", width, height);
  return SpzFrameWritePlan{std::move(data),
                           [header = encodeHeader(*this)](std::ostream& os, const std::vector<SpzFrameLocation>& frames) {
                             std::ostringstream payload;
                             payload << header;
                             writeFrameTable(payload, frames);
                             writeExtensionRecord(os, SpzExtensionType::SPZ_ADOBE_environment_map, payload.str());
                           }};
}

std::optional<SpzFrameReadPlan> parseEnvironmentMapAdobe(std::istream& is) {
  SpzExtensionEnvironmentMapAdobe header;
  uint32_t metadataLength{};
  bool ok = readExact(is, header.width) && readExact(is, header.height);
  for (auto* v : {&header.origin, &header.up, &header.e1, &header.e2})
    ok = ok && is.read(reinterpret_cast<char*>(v->data()), v->size() * sizeof(float));
  ok = ok && readExact(is, metadataLength);
  // The payload stream is an in-memory buffer, so in_avail() is the rest of it.
  ok = ok && metadataLength <= is.rdbuf()->in_avail();
  if (ok) {
    header.metadata.resize(metadataLength);
    ok = static_cast<bool>(is.read(header.metadata.data(), metadataLength));
  }
  std::optional<std::vector<SpzFrameLocation>> frames = ok ? readFrameTable(is) : std::nullopt;
  if (!frames) {
    SpzLog("[SPZ ERROR] EnvironmentMapAdobe: truncated record");
    return std::nullopt;
  }
  if (header.width == 0 || header.height == 0) {
    SpzLog("[SPZ ERROR] EnvironmentMapAdobe: width and height must be at least 1");
    return std::nullopt;
  }
  SpzLog("[SPZ] Found extension: EnvironmentMapAdobe (%ux%u)", header.width, header.height);
  const uint32_t width = header.width;
  const uint32_t height = header.height;
  return SpzFrameReadPlan{
    SpzExtensionType::SPZ_ADOBE_environment_map, std::move(*frames), header.uniquenessKey(),
    [width, height](int32_t) -> std::optional<uint64_t> { return pixelBytes(width, height); },
    [header = std::move(header)](int32_t, std::vector<uint8_t>&& data) -> SpzExtensionBasePtr {
      auto rec = std::make_shared<SpzExtensionEnvironmentMapAdobe>(header);
      rec->pixels.resize(data.size() / sizeof(Half));
      for (size_t i = 0; i < rec->pixels.size(); ++i) {
        Half half;
        std::memcpy(&half, &data[i * sizeof(Half)], sizeof(Half));
        rec->pixels[i] = halfToFloat(half);
      }
      if (!rec->isValid()) return nullptr;
      return rec;
    }};
}

SpzExtensionType SpzExtensionEnvironmentMapAdobe::type() {
  return SpzExtensionType::SPZ_ADOBE_environment_map;
}

SpzExtensionBase* SpzExtensionEnvironmentMapAdobe::copyAsRawData() const {
  return new SpzExtensionEnvironmentMapAdobe(*this);
}

std::optional<SpzExtensionBasePtr> SpzExtensionEnvironmentMapAdobe::tryReadFromPly(
    std::istream& in, const std::vector<PlyExtraElement>& elements) const {
  if (elements.size() != kRequiredPlyElementNames.size() || elements[0].name != kPlyEnvironment ||
      elements[0].count != 1 || elements[0].bytesPerElement != kFixedBytes - sizeof(uint32_t) ||
      elements[1].name != kPlyMap || elements[1].bytesPerElement != kChannels * sizeof(float) ||
      elements[2].name != kPlyMetadata || elements[2].bytesPerElement != 1) {
    SpzLog("[SPZ ERROR] Environment map: PLY elements don't match the expected layout");
    return std::nullopt;
  }
  auto rec = std::make_shared<SpzExtensionEnvironmentMapAdobe>();
  bool ok = readExact(in, rec->width) && readExact(in, rec->height);
  for (auto* v : {&rec->origin, &rec->up, &rec->e1, &rec->e2})
    ok = ok && in.read(reinterpret_cast<char*>(v->data()), 3 * sizeof(float));
  if (!ok || static_cast<uint64_t>(elements[1].count) != static_cast<uint64_t>(rec->width) * rec->height) {
    SpzLog("[SPZ ERROR] Environment map: PLY header doesn't match %ux%u", rec->width, rec->height);
    return std::nullopt;
  }
  rec->pixels.resize(pixelFloats(rec->width, rec->height));
  rec->metadata.resize(static_cast<size_t>(elements[2].count));
  if (!in.read(reinterpret_cast<char*>(rec->pixels.data()), rec->pixels.size() * sizeof(float)) ||
      !in.read(rec->metadata.data(), static_cast<std::streamsize>(rec->metadata.size()))) {
    SpzLog("[SPZ ERROR] Environment map: truncated PLY data");
    return std::nullopt;
  }
  if (!rec->isValid())
    return std::nullopt;
  SpzLog("[SPZ] Loaded environment map: %ux%u", rec->width, rec->height);
  return std::optional{std::move(rec)};
}

void SpzExtensionEnvironmentMapAdobe::writePlyHeader(std::ostream& out) const {
  out << "element " << kPlyEnvironment << " 1\n";
  out << "property uint environment_width\n";
  out << "property uint environment_height\n";
  for (const char* v : {"origin", "up", "e1", "e2"}) {
    for (const char* axis : {"x", "y", "z"})
      out << "property float environment_" << v << "_" << axis << "\n";
  }
  out << "element " << kPlyMap << " " << static_cast<uint64_t>(width) * height << "\n";
  for (const char* channel : {"r", "g", "b"})
    out << "property float " << kPlyMap << "_" << channel << "\n";
  out << "element " << kPlyMetadata << " " << metadata.size() << "\n";
  out << "property uchar " << kPlyMetadata << "\n";
}

void SpzExtensionEnvironmentMapAdobe::writePlyData(std::ostream& out) const {
  out.write(reinterpret_cast<const char*>(&width), sizeof(width));
  out.write(reinterpret_cast<const char*>(&height), sizeof(height));
  for (const auto* v : {&origin, &up, &e1, &e2})
    out.write(reinterpret_cast<const char*>(v->data()), 3 * sizeof(float));
  out.write(reinterpret_cast<const char*>(pixels.data()), pixels.size() * sizeof(float));
  out.write(metadata.data(), static_cast<std::streamsize>(metadata.size()));
}

}  // namespace spz
