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
#include "load-spz.h"

#include <zstd.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace spz {
namespace {
constexpr int kZstdLevel = 3;
// Same bound as custom attributes: no zstd frame inflates more than 32768x.
constexpr uint64_t kMaxZstdExpansion = 32768;
constexpr float kMaxHalf = 65504.0f;
// width, height, origin, up, e1, e2, metadataLength
constexpr size_t kFixedBytes = 2 * sizeof(uint32_t) + 12 * sizeof(float) + sizeof(uint32_t);

template <class T>
bool readExact(std::istream& is, T& out) {
  return static_cast<bool>(is.read(reinterpret_cast<char*>(&out), sizeof(T)));
}

template <class T>
void append(std::string& out, const T* data, size_t count) {
  out.append(reinterpret_cast<const char*>(data), count * sizeof(T));
}

bool isFinite(const std::array<float, 3>& v) {
  return std::all_of(v.begin(), v.end(), [](float x) { return std::isfinite(x); });
}

uint64_t pixelFloats(uint32_t width, uint32_t height) {
  return static_cast<uint64_t>(width) * height * 3;
}
}  // namespace

const std::unordered_set<std::string> SpzExtensionEnvironmentMapAdobe::kRequiredPlyElementNames = {
  "environment",
  "environment_map",
  "environment_metadata",
};

SpzExtensionEnvironmentMapAdobe::SpzExtensionEnvironmentMapAdobe()
    : SpzExtensionBase(SpzExtensionType::SPZ_ADOBE_environment_map) {}

bool SpzExtensionEnvironmentMapAdobe::validate() const {
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
  else if (kFixedBytes + metadata.size() + ZSTD_compressBound(pixels.size() * sizeof(Half)) >
           std::numeric_limits<uint32_t>::max())
    problem = "data exceeds the 4 GiB record limit";
  if (problem) {
    SpzLog("[SPZ ERROR] Environment map: %s", problem);
    return false;
  }
  return true;
}

std::string SpzExtensionEnvironmentMapAdobe::encode() const {
  const uint32_t metadataLength = static_cast<uint32_t>(metadata.size());
  std::string out;
  append(out, &width, 1);
  append(out, &height, 1);
  for (const auto* v : {&origin, &up, &e1, &e2})
    append(out, v->data(), 3);
  append(out, &metadataLength, 1);
  out += metadata;

  std::vector<Half> halves(pixels.size());
  std::transform(pixels.begin(), pixels.end(), halves.begin(), floatToHalf);
  const size_t srcBytes = halves.size() * sizeof(Half);
  const size_t offset = out.size();
  out.resize(offset + ZSTD_compressBound(srcBytes));
  size_t compressed = ZSTD_compress(&out[offset], out.size() - offset, halves.data(), srcBytes, kZstdLevel);
  if (ZSTD_isError(compressed)) {
    // Writing no pixels makes readers reject this record instead of misreading it.
    SpzLog("[SPZ ERROR] Environment map: zstd compression failed");
    compressed = 0;
  }
  out.resize(offset + compressed);
  return out;
}

uint32_t SpzExtensionEnvironmentMapAdobe::payloadBytes() const {
  return static_cast<uint32_t>(encode().size());
}

void SpzExtensionEnvironmentMapAdobe::write(std::ostream& os) const {
  const std::string payload = encode();
  const uint32_t t = static_cast<uint32_t>(extensionType);
  const uint32_t len = static_cast<uint32_t>(payload.size());  // validate() keeps this under 4 GiB
  SpzLog("[SPZ] Writing extension: EnvironmentMapAdobe (%ux%u, %u bytes)", width, height, len);
  os.write(reinterpret_cast<const char*>(&t), sizeof(t));
  os.write(reinterpret_cast<const char*>(&len), sizeof(len));
  os.write(payload.data(), static_cast<std::streamsize>(payload.size()));
}

std::optional<SpzExtensionBasePtr> SpzExtensionEnvironmentMapAdobe::read(std::istream& is) {
  auto rec = std::make_shared<SpzExtensionEnvironmentMapAdobe>();
  uint32_t metadataLength{};
  bool ok = readExact(is, rec->width) && readExact(is, rec->height);
  for (auto* v : {&rec->origin, &rec->up, &rec->e1, &rec->e2})
    ok = ok && is.read(reinterpret_cast<char*>(v->data()), 3 * sizeof(float));
  ok = ok && readExact(is, metadataLength);
  // The payload stream is an in-memory buffer, so in_avail() is the rest of it.
  ok = ok && metadataLength < is.rdbuf()->in_avail();
  if (ok) {
    rec->metadata.resize(metadataLength);
    ok = static_cast<bool>(is.read(rec->metadata.data(), metadataLength));
  }
  if (!ok) {
    SpzLog("[SPZ ERROR] EnvironmentMapAdobe: truncated payload");
    return std::nullopt;
  }

  const std::streamsize compressedSize = is.rdbuf()->in_avail();
  const uint64_t expected = pixelFloats(rec->width, rec->height) * sizeof(Half);
  if (expected == 0 || expected > std::numeric_limits<size_t>::max() ||
      expected > static_cast<uint64_t>(compressedSize) * kMaxZstdExpansion) {
    SpzLog("[SPZ ERROR] EnvironmentMapAdobe: pixels don't cover %ux%u", rec->width, rec->height);
    return std::nullopt;
  }
  std::string compressed(static_cast<size_t>(compressedSize), '\0');
  is.read(compressed.data(), compressedSize);
  std::vector<Half> halves(static_cast<size_t>(expected / sizeof(Half)));
  const size_t decompressed =
      ZSTD_decompress(halves.data(), static_cast<size_t>(expected), compressed.data(), compressed.size());
  if (ZSTD_isError(decompressed) || decompressed != expected) {
    SpzLog("[SPZ ERROR] EnvironmentMapAdobe: pixels don't cover %ux%u", rec->width, rec->height);
    return std::nullopt;
  }
  rec->pixels.resize(halves.size());
  std::transform(halves.begin(), halves.end(), rec->pixels.begin(), halfToFloat);
  if (!rec->validate())
    return std::nullopt;

  SpzLog("[SPZ] Found extension: EnvironmentMapAdobe (%ux%u)", rec->width, rec->height);
  return std::optional{std::move(rec)};
}

SpzExtensionType SpzExtensionEnvironmentMapAdobe::type() {
  return SpzExtensionType::SPZ_ADOBE_environment_map;
}

SpzExtensionBase* SpzExtensionEnvironmentMapAdobe::copyAsRawData() const {
  return new SpzExtensionEnvironmentMapAdobe(*this);
}

std::optional<SpzExtensionBasePtr> SpzExtensionEnvironmentMapAdobe::tryReadFromPly(
    std::istream& in, const std::vector<PlyExtraElement>& elements) const {
  if (elements.size() != 3 || elements[0].name != "environment" || elements[0].count != 1 ||
      elements[0].bytesPerElement != kFixedBytes - sizeof(uint32_t) ||
      elements[1].name != "environment_map" || elements[1].bytesPerElement != 3 * sizeof(float) ||
      elements[2].name != "environment_metadata" || elements[2].bytesPerElement != 1) {
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
  if (!rec->validate())
    return std::nullopt;
  SpzLog("[SPZ] Loaded environment map: %ux%u", rec->width, rec->height);
  return std::optional{std::move(rec)};
}

void SpzExtensionEnvironmentMapAdobe::writePlyHeader(std::ostream& out) const {
  out << "element environment 1\n";
  out << "property uint environment_width\n";
  out << "property uint environment_height\n";
  for (const char* v : {"origin", "up", "e1", "e2"}) {
    for (const char* axis : {"x", "y", "z"})
      out << "property float environment_" << v << "_" << axis << "\n";
  }
  out << "element environment_map " << static_cast<uint64_t>(width) * height << "\n";
  out << "property float environment_map_r\n";
  out << "property float environment_map_g\n";
  out << "property float environment_map_b\n";
  out << "element environment_metadata " << metadata.size() << "\n";
  out << "property uchar environment_metadata\n";
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
