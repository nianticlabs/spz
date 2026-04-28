# SPZ Extensions

Extensions allow vendor-specific or application-specific data to be stored in SPZ files alongside the core Gaussian splat data. Multiple extensions from different vendors can coexist in the same file; unknown extension types are skipped during parsing so readers only need to understand the extensions they care about.

## Extension stream format

Extensions are stored at the end of the packed SPZ stream (after Gaussian data). The format is **length‑delimited per record**:

```
[ u32 type ][ u32 byteLength ][ payload... ]
```

- **type** — 4-byte extension type ID (e.g. `0xADBE0002` for Adobe safe orbit camera).
- **byteLength** — 4-byte length in bytes of the following payload.
- **payload** — exactly `byteLength` bytes of extension-specific data.

Records repeat back-to-back until end-of-stream. There is no separate record count or terminator; EOF after the last payload ends the extension block.

### Why length-delimited?

- **Unknown types can be skipped**: a reader that doesn’t recognize `type` can skip `byteLength` bytes and continue to the next record.
- **Multiple vendors**: different vendors can use their own type IDs and payloads in the same file without conflicting.

When the library encounters an **unknown extension type** (a type it does not implement), it skips that record and continues: it logs a warning including the type ID and payload size (e.g. `[SPZ WARNING] Skipping unknown extension type: 0x12340001 (24 bytes)`), advances the stream past the payload by `byteLength` bytes, and parses the next record. The load succeeds; only extensions the library knows are added to the result. This allows files to carry extensions from multiple vendors while each reader uses only the ones it supports.

### Extension type ID allocation (avoiding collisions)

The 32-bit extension type is split into two 16-bit fields so multiple vendors can allocate IDs without coordination:

| Bits    | Field         | Range (hex) | Meaning |
|--------|----------------|------------|---------|
| 31–16  | **Vendor ID**  | `0x0001`–`0xFFFF` | Unique per vendor; pick a value and use it for all your extensions. |
| 15–0   | **Extension ID** | `0x0001`–`0xFFFF` | Unique within that vendor; you can use any numbering scheme. |

- **Layout:** `type_u32 == (VENDOR_ID << 16) | EXTENSION_ID`.  
- **Reserved:** `0x00000000` is not used as an extension type (no vendor 0).  
- **Vendor ID:** Choose a 16-bit value that is unique to your organization or product (e.g. Adobe uses `0xADBE`). To avoid collisions, avoid well-known prefixes (e.g. `0xADBE`, `0x4E41` for Niantic) and document or register your vendor ID if you publish extensions.  
- **Extension ID:** Within your vendor space, use the low 16 bits for each extension (e.g. `0x0001`, `0x0002`, …). No global registry; uniqueness per vendor is enough.

Example: Adobe vendor ID `0xADBE`, extension index 2 → `0xADBE0002u` (`SPZ_ADOBE_safe_orbit_camera`).

## Builds without extension support

If an SPZ file that contains extensions is loaded by a library build that was **not** built with extension support (`SPZ_BUILD_EXTENSIONS` is OFF, the default):

- **The load still succeeds.** Core Gaussian splat data (positions, colors, SH, etc.) is read and returned as usual.
- **Extension data is ignored.** The header’s extension flag is read, but the extension block is not parsed or skipped. A warning is logged:  
  `[SPZ WARNING] deserializePackedGaussians: the stream has extensions but extensions are unsupported in the current build of SPZ`
- **Returned clouds have no extensions.** `PackedGaussians` and `GaussianCloud` in a no-extension build do not have an `extensions` member; the loaded result simply omits any extension data.

So builds without extensions can safely load files that contain extensions; they get the core data and drop the rest. To preserve or use extension data, build with `-DSPZ_BUILD_EXTENSIONS=ON` and use a build that includes the extension sources.

## Using extensions

### C++

1. **Build with extensions enabled**  
   CMake: `-DSPZ_BUILD_EXTENSIONS=ON`.

2. **Load/Save**  
   Extensions are part of `PackedGaussians` and `GaussianCloud`. When you load (e.g. `loadSpz`, `loadSpzPacked`) or save (`saveSpz`, `serializePackedGaussians`), the extension list is read/written automatically.

3. **Access a known extension**  
   Use `findExtensionByType<T>()`. Include the header for the concrete extension type you use:

   ```cpp
   #include "splat-extensions.h"
   #include "safe-orbit-camera-adobe.h"

   auto ext = spz::findExtensionByType<spz::SpzExtensionSafeOrbitCameraAdobe>(cloud.extensions);
   if (ext) {
     float minElev = ext->safeOrbitElevationMin;
     float maxElev = ext->safeOrbitElevationMax;
     float minRadius = ext->safeOrbitRadiusMin;
   }
   ```

4. **Add an extension before saving**  
   Push a shared pointer onto the extensions vector. Include the extension’s header:

   ```cpp
   #include "splat-extensions.h"
   #include "safe-orbit-camera-adobe.h"

   auto ext = std::make_shared<spz::SpzExtensionSafeOrbitCameraAdobe>();
   ext->safeOrbitElevationMin = -0.5f;
   ext->safeOrbitElevationMax = 0.5f;
   ext->safeOrbitRadiusMin = 1.0f;
   packed.extensions.push_back(ext);
   ```

### Python

With extension support built in, `GaussianCloud.extensions` is a list of extension objects. Check `extension_type` and cast or use type-specific classes (e.g. `SpzExtensionSafeOrbitCameraAdobe`) as exposed by the bindings.

### JavaScript / TypeScript (WASM)

The WASM build exposes extension types and the cloud’s `extensions` array. Use the generated TypeScript types and the same pattern: iterate or filter by `extensionType`, then use the concrete extension interface.

### PLY extension element names

When loading PLY files, the loader treats non-vertex elements as either “known” (handled by an extension) or unknown. Known elements are consumed by `readExtensionsFromPly`; unknown ones are reported and skipped. The set of known names is derived from the **PLY extension registry** in `extensions/cc/splat-extensions.cc`: each extension that supports PLY declares a static member (e.g. `SpzExtensionSafeOrbitCameraAdobe::kRequiredPlyElementNames`) and is registered with an exemplar in `getPlyExtensionRegistry()`. The union of all registered extension element names is used to answer **`isKnownPlyExtensionElement(const std::string& elementName)`**. To add PLY support for a new extension, implement the three abstract PLY methods on your extension struct and register it in the registry (see “Adding a new extension” below).

## Adding a new extension

To add a new extension type in the C++ codebase:

1. **Define the type ID**  
   In `extensions/cc/splat-extensions.h`, add a value to `enum class SpzExtensionType` using your vendor prefix, e.g.:

   ```cpp
   enum class SpzExtensionType : uint32_t {
     SPZ_ADOBE_safe_orbit_camera = 0xADBE0002u,
     SPZ_MYVENDOR_my_extension    = 0xXXXX0001u,  // your vendor ID and ID
   };
   ```

2. **Define the extension struct**  
   Declare a struct that inherits `SpzExtensionBase`, with your payload fields and the required overrides. You can add it in `extensions/cc/splat-extensions.h` or in a separate header (e.g. `my-extension.h`) that `#include "splat-extensions.h"`; the built-in Adobe safe orbit extension lives in `safe-orbit-camera-adobe.h` / `safe-orbit-camera-adobe.cc` as a reference. The struct must have:

   - Constructor.
   - `uint32_t payloadBytes() const override` — return the payload size in bytes.
   - `void write(std::ostream& os) const override` — emit type (4 bytes), byteLength (4 bytes), then payload; use `payloadBytes()` for the length (see step 3).
   - `SpzExtensionBase* copyAsRawData() const override` — allocate a copy (e.g. `return new MyExtension(*this);`).
   - `static std::optional<SpzExtensionBasePtr> read(std::istream& is)` — read payload from `is` (caller provides a stream of exactly your payload bytes).
   - `static SpzExtensionType type()` — return your enum value.
   - The three PLY methods (required because they are pure virtual on the base): `tryReadFromPly(std::istream& in, const std::unordered_set<std::string>& elementNames) const override`, `writePlyHeader(std::ostream& out) const override`, `writePlyData(std::ostream& out) const override`. If your extension does not support PLY, implement them to return `std::nullopt` / no-op as appropriate; if it does, see step 6.

3. **Implement write() and payloadBytes()**  
   Implement `payloadBytes()` to return your payload size. Implement `write()` to emit type, byteLength (use `payloadBytes()`), then payload. Example:

   ```cpp
   uint32_t MyExtension::payloadBytes() const {
     return sizeof(field1) + sizeof(field2);  // your payload size
   }

   void MyExtension::write(std::ostream& os) const {
     const uint32_t t = static_cast<uint32_t>(extensionType);
     const uint32_t len = payloadBytes();
     os.write(reinterpret_cast<const char*>(&t), sizeof(t));
     os.write(reinterpret_cast<const char*>(&len), sizeof(len));
     os.write(reinterpret_cast<const char*>(&field1), sizeof(field1));
     // ... rest of payload
   }
   ```

4. **Implement read()**  
   Your `read(std::istream& is)` receives a stream positioned at the **start of your payload** (the caller has already read type and byteLength and verified the length). Read exactly your payload and return an `std::optional<SpzExtensionBasePtr>`. On parse failure, return `std::nullopt` (if the stream is a temporary buffer, e.g. from an extension payload, no rollback is needed).

5. **Register in the parser**  
   In `extensions/cc/splat-extensions.cc`, ensure your extension’s header is included (e.g. `#include "my-extension.h"`). In `tryParseExtension`, add a `case` for your `SpzExtensionType`: read `byteLength` bytes into a buffer, open an `std::istringstream` on it, and call `MyExtension::read(iss)`. Push the result into `out` if valid. Return `true` to continue. Unknown types are already handled: the common code skips `byteLength` bytes and continues.

6. **Optional: PLY round-trip**  
   If you want your extension to load/save from PLY as well:

   - **Required PLY element names** — In your extension struct, add a static member, e.g. `static const std::unordered_set<std::string> kRequiredPlyElementNames;`, and define it in the .cc with the PLY extra element names your extension uses.
   - **Implement the three PLY methods** — `tryReadFromPly(in, elementNames)`: if `elementNames` contains all of `kRequiredPlyElementNames`, read your bytes from `in` and return a new instance (or `std::nullopt` on failure). `writePlyHeader(out)`: write the PLY `element` / `property` lines for your data. `writePlyData(out)`: write the binary data for your extension.
   - **Register in the PLY registry** — In `extensions/cc/splat-extensions.cc`, inside `getPlyExtensionRegistry()`, add an entry: `{&MyExtension::kRequiredPlyElementNames, std::make_shared<MyExtension>()}`. Order in the registry is the order used when reading/writing PLY. You do **not** edit `readExtensionsFromPly`, `writeExtensionsToPlyHeader`, or `writeExtensionsToPlyData`; they loop over the registry and call your virtual methods.

7. **Bindings**  
   Expose the new type and struct in Python (`extensions/python/splat-extensions.cc`) and Emscripten (`extensions/emscripten/splat-extensions.cc`, plus `.d.ts.in` if needed) so callers can construct and read your extension from the cloud’s `extensions` list. Each binding file must `#include` your extension’s header (e.g. `#include "safe-orbit-camera-adobe.h"`) to get the full type definition.

## Built-in extensions

- **SPZ_ADOBE_safe_orbit_camera** (`0xADBE0002`) — Camera orbit limits (elevation min/max, radius min) for restricting the view. Implemented in `extensions/cc/safe-orbit-camera-adobe.h` and `safe-orbit-camera-adobe.cc`. See the main [README](../README.md) for attributes and defaults.
