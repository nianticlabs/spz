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
#include <functional>
#include <iostream>
#include <istream>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <unordered_set>
#include <vector>

// C-compatible extension node structure (used for C interop)
typedef struct SpzExtensionNode {
  uint32_t type;
  void* data;
  struct SpzExtensionNode* next;
} SpzExtensionNode;

namespace spz {
// Forward declaration to avoid circular dependency
// The full definition is in load-spz.h, which is included in splat-extensions.cc
struct PlyExtraElement;
struct CoordinateConverter;
enum class CoordinateSystem : uint32_t;
enum class SpzExtensionType : uint32_t {
  SPZ_ADOBE_safe_orbit_camera = 0xADBE0002u,
  SPZ_ADOBE_coordinate_system = 0xADBE0003u,
  SPZ_NIANTIC_georeference = 0x4E530001u,
  SPZ_NIANTIC_custom_attribute = 0x4E530002u,
  SPZ_ADOBE_environment_map = 0xADBE0005u,
};

struct SpzExtensionBase {
  SpzExtensionType extensionType;
  explicit SpzExtensionBase(SpzExtensionType t) : extensionType(t) {
  }
  virtual ~SpzExtensionBase() = default;
  virtual uint32_t payloadBytes() const = 0;
  virtual void write(std::ostream& os) const = 0;
  // False when write() would emit a record that read() rejects; saveSpz then refuses to save.
  virtual bool isValid() const { return true; }
  virtual SpzExtensionBase* copyAsRawData() const = 0;

  // PLY extension I/O: read from PLY, write header lines, write data. `elements` is the run of
  // consecutive PLY elements named in the extension's kRequiredPlyElementNames, in file order;
  // `in` is at the start of their data.
  virtual std::optional<std::shared_ptr<SpzExtensionBase>> tryReadFromPly(
      std::istream& in, const std::vector<PlyExtraElement>& elements) const = 0;
  virtual void writePlyHeader(std::ostream& out) const = 0;
  virtual void writePlyData(std::ostream& out) const = 0;
};

using SpzExtensionBasePtr = std::shared_ptr<SpzExtensionBase>;

// Stream helpers shared by extension readers.
namespace detail {
template <class T>
bool readExact(std::istream& is, T& out) {
  return static_cast<bool>(is.read(reinterpret_cast<char*>(&out), sizeof(T)));
}

// Returns the number of bytes remaining between the current get position and the end of `is`,
// or -1 if the stream does not support positioning. Restores the original position on success.
std::streamoff remainingBytes(std::istream& is);
}  // namespace detail

// ---- Optional capabilities ----
// An extension opts into a capability by also inheriting its interface; the library detects it
// with dynamic_cast. Capabilities are separate from SpzExtensionBase so extensions built against an
// older SpzExtensionBase keep working unchanged.

// ---- Extension frames: storing bulk data outside the extension zone ----
// An extension frame is a zstd skippable frame (magic 0x184D2A53) whose payload is exactly one zstd
// frame holding part of an extension's data, appended after the last core stream's own zstd frame
// in an NGSP (v4) file (see "Extension frames" in README.md).

// Location of one extension frame. frameOffset is relative to the start of the last core stream's
// data and points at the magic number. compressedSize is the payload size (excluding the 8-byte
// skippable header); uncompressedSize is the payload's decompressed size.
struct SpzFrameLocation {
  uint64_t frameOffset = 0;
  uint64_t compressedSize = 0;
  uint64_t uncompressedSize = 0;
};

// How to save a framed record: its bulk data, and how to write the record once the library has
// placed that data in extension frames. Returned by SpzFramedExtension::planFrameWrites.
struct SpzFrameWritePlan {
  // The record's bulk data. The library splits it across extension frames and compresses them.
  std::vector<uint8_t> data;
  // Writes the whole record (type, byteLength, payload). Precondition: frames holds one location
  // per extension frame of data, in order.
  std::function<void(std::ostream& os, const std::vector<SpzFrameLocation>& frames)> writeRecord;
};

// Capability: the record keeps its bulk data in extension frames instead of in the plaintext
// extension zone. Readers that don't know the extension skip the extension frames without
// decompressing them, so files stay loadable by every v4 reader.
struct SpzFramedExtension {
  virtual ~SpzFramedExtension() = default;
  // Returns the save plan for a cloud of numPoints points, or std::nullopt (after logging why) if
  // the record is invalid and the save must fail.
  virtual std::optional<SpzFrameWritePlan> planFrameWrites(int32_t numPoints) const = 0;
};


// Capability: no two records of the same type may share a key (e.g. an attribute's name). saveSpz
// fails if two do, and loading drops every record that shares a key, since none can be told apart.
struct SpzKeyedExtension {
  virtual ~SpzKeyedExtension() = default;
  virtual std::string uniquenessKey() const = 0;
};

// Capability: the record holds values expressed in the cloud's coordinate system, which must
// change along with the cloud's coordinates.
struct SpzCoordinateDependentExtension {
  virtual ~SpzCoordinateDependentExtension() = default;
  // Returns a converted copy, or nullptr if c leaves the record unchanged. Records are shared
  // between clouds, so implementations must never modify *this.
  virtual SpzExtensionBasePtr withConvertedCoordinates(const CoordinateConverter& c) const = 0;
};

// Built-in extension types are declared in their own headers, e.g. safe-orbit-camera-adobe.h
// Reads the records that need no frame data; framed records are dropped with an error, since only
// NGSP (v4) files can carry their data.
void readAllExtensions(std::istream& is, std::vector<SpzExtensionBasePtr>& out);

// Writes every record that needs no frame data. Framed records (SpzFramedExtension) are skipped
// with an error; saveSpz writes them through their SpzFrameWritePlan.
void writeAllExtensions(const std::vector<SpzExtensionBasePtr>& list, std::ostream& os);

void readExtensionsFromPly(std::istream& in, const std::vector<spz::PlyExtraElement>& extraElements, std::vector<SpzExtensionBasePtr>& extensions);

void writeExtensionsToPlyHeader(const std::vector<SpzExtensionBasePtr>& extensions, std::ostream& out);

void writeExtensionsToPlyData(const std::vector<SpzExtensionBasePtr>& extensions, std::ostream& out);

bool isKnownPlyExtensionElement(const std::string& elementName);

// Returns the coordinate system in which the Gaussian data is stored, as recorded by any
// behavior-modifying extension. Falls back to RUB if no such extension is present.
CoordinateSystem getPackedCoordinateSystem(const std::vector<SpzExtensionBasePtr>& extensions);

// False if any extension reports isValid() == false; the extension logs the reason.
bool allExtensionsValid(const std::vector<SpzExtensionBasePtr>& extensions);

inline SpzExtensionNode* copyExtensions(const std::vector<SpzExtensionBasePtr> &extensions) {
  SpzExtensionNode* head = nullptr;
  SpzExtensionNode* tail = nullptr;

  for (const auto& ext : extensions) {
    if (!ext) continue;  // Skip null extensions

    SpzExtensionNode* node = new SpzExtensionNode{static_cast<uint32_t>(ext->extensionType), ext->copyAsRawData(), nullptr};

    if (!head) {
      head = node;
      tail = node;
    } else {
      tail->next = node;
      tail = node;
    }
  }

  return head;
}

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
