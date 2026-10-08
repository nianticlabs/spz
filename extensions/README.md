# SPZ Extensions

Extensions allow vendor-specific or application-specific data to be stored in SPZ files alongside the core Gaussian splat data. Multiple extensions from different vendors can coexist in the same file; unknown extension types are skipped during parsing so readers only need to understand the extensions they care about.

## Extension record format

Extensions use a **length-delimited per-record** format regardless of where they are stored:

```
[ u32 type ][ u32 byteLength ][ payload... ]
```

- **type** — 4-byte extension type ID (e.g. `0xADBE0002` for Adobe safe orbit camera).
- **byteLength** — 4-byte length in bytes of the following payload.
- **payload** — exactly `byteLength` bytes of extension-specific data.

Records repeat back-to-back. Extensions occupy the plaintext header zone between the fixed 32-byte `NgspFileHeader` and the Table of Contents. The block spans bytes `[32, tocByteOffset)`, where `tocByteOffset` is read directly from the header — no need to decompress anything to locate or parse extensions.

### Why length-delimited?

- **Unknown types can be skipped**: a reader that doesn’t recognize `type` can skip `byteLength` bytes and continue to the next record.
- **Multiple vendors**: different vendors can use their own type IDs and payloads in the same file without conflicting.

When the library encounters an **unknown extension type** (a type it does not implement), it skips that record and continues: it logs a warning including the type ID and payload size (e.g. `[SPZ WARNING] Unknown extension type 0x12340001 (24 bytes) was skipped — loaded data may be incorrect`), advances the stream past the payload by `byteLength` bytes, and parses the next record. The load succeeds; only extensions the library knows are added to the result. This allows files to carry extensions from multiple vendors while each reader uses only the ones it supports.

### Extension type ID allocation (avoiding collisions)

The 32-bit extension type is split into two 16-bit fields so multiple vendors can allocate IDs without coordination:

| Bits    | Field         | Range (hex) | Meaning |
|--------|----------------|------------|---------|
| 31–16  | **Vendor ID**  | `0x0001`–`0xFFFF` | Unique per vendor; pick a value and use it for all your extensions. |
| 15–0   | **Extension ID** | `0x0001`–`0xFFFF` | Unique within that vendor; you can use any numbering scheme. |

- **Layout:** `type_u32 == (VENDOR_ID << 16) | EXTENSION_ID`.  
- **Reserved:** `0x00000000` is not used as an extension type (no vendor 0).  
- **Vendor ID:** Choose a 16-bit value that is unique to your organization or product (e.g. Adobe uses `0xADBE`). To avoid collisions, avoid the allocated vendor IDs below and document or register your vendor ID if you publish extensions.  
- **Extension ID:** Within your vendor space, use the low 16 bits for each extension (e.g. `0x0001`, `0x0002`, …). No global registry; uniqueness per vendor is enough.

Allocated vendor IDs:

| Vendor ID | Vendor |
|-----------|--------|
| `0xADBE`  | Adobe |
| `0x4E53`  | Niantic (`'NS'`) |

Example: Adobe vendor ID `0xADBE`, extension index 2 → `0xADBE0002u` (`SPZ_ADOBE_safe_orbit_camera`).

## Builds without extension support

If an SPZ file that contains extensions is loaded by a library build that was **not** built with extension support (`SPZ_BUILD_EXTENSIONS` is OFF, the default):

- **The load still succeeds.** Core Gaussian splat data (positions, colors, SH, etc.) is read and returned as usual.
- **Extension data is ignored.** The header’s extension flag is read, but the extension block is not parsed. A warning is logged. The exact message depends on the file format:
  - NGSP v4: `[SPZ WARNING] loadSpzPacked: file has extensions but extension support is disabled — skipped extensions may affect how data was packed or will be unpacked; build with SPZ_BUILD_EXTENSIONS to ensure correct results`
  - Legacy gzip: `[SPZ WARNING] deserializePackedGaussians: stream has extensions but extension support is disabled — skipped extensions may affect how data was packed or will be unpacked; build with SPZ_BUILD_EXTENSIONS to ensure correct results`
- **Returned clouds have no extensions.** `PackedGaussians` and `GaussianCloud` in a no-extension build do not have an `extensions` member; the loaded result simply omits any extension data.
- **Data may be incorrect.** Some extensions (such as `SPZ_ADOBE_coordinate_system`) alter how data is packed or unpacked. If a file was written with such an extension, a build without extension support will silently misinterpret the data — for example, applying the wrong coordinate system conversion. The warning above is the only signal this has occurred.

To preserve or use extension data, and to ensure correct data interpretation, build with `-DSPZ_BUILD_EXTENSIONS=ON`.

## Using extensions

### C++

1. **Build with extensions enabled**  
   CMake: `-DSPZ_BUILD_EXTENSIONS=ON`.

2. **Load/Save**  
   Extensions are part of `PackedGaussians` and `GaussianCloud`. When you call `loadSpz` / `loadSpzPacked` or `saveSpz`, the extension list is read/written automatically. For version 4 and above, `saveSpz` writes extensions into the plaintext header zone.

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
   In `extensions/cc/splat-extensions.cc`, ensure your extension’s header is included (e.g. `#include "my-extension.h"`). In `tryParseExtension`, add a `case` for your `SpzExtensionType`: read `byteLength` bytes into a buffer, open an `std::istringstream` on it, and call `MyExtension::read(iss)`. Push the result into `out` if valid. Return `true` to continue. Unknown types are already handled: the common code skips `byteLength` bytes and logs a warning.

   If your extension affects the packed coordinate system, also add a case in `getPackedCoordinateSystem` (in the same file) that calls your extension’s resolve logic. This ensures `packGaussians` and `unpackGaussians` apply the correct coordinate conversion.

6. **Optional: PLY round-trip**  
   If you want your extension to load/save from PLY as well:

   - **Required PLY element names** — In your extension struct, add a static member, e.g. `static const std::unordered_set<std::string> kRequiredPlyElementNames;`, and define it in the .cc with the PLY extra element names your extension uses.
   - **Implement the three PLY methods** — `tryReadFromPly(in, elementNames)`: if `elementNames` contains all of `kRequiredPlyElementNames`, read your bytes from `in` and return a new instance (or `std::nullopt` on failure). `writePlyHeader(out)`: write the PLY `element` / `property` lines for your data. `writePlyData(out)`: write the binary data for your extension.
   - **Register in the PLY registry** — In `extensions/cc/splat-extensions.cc`, inside `getPlyExtensionRegistry()`, add an entry: `{&MyExtension::kRequiredPlyElementNames, std::make_shared<MyExtension>()}`. Order in the registry is the order used when reading/writing PLY. You do **not** edit `readExtensionsFromPly`, `writeExtensionsToPlyHeader`, or `writeExtensionsToPlyData`; they loop over the registry and call your virtual methods.

7. **Bindings**  
   Expose the new type and struct in Python (`extensions/python/splat-extensions.cc`) and Emscripten (`extensions/emscripten/splat-extensions.cc`, plus `.d.ts.in` if needed) so callers can construct and read your extension from the cloud’s `extensions` list. Each binding file must `#include` your extension’s header (e.g. `#include "safe-orbit-camera-adobe.h"`) to get the full type definition.

## Built-in extensions

- **SPZ_ADOBE_safe_orbit_camera** (`0xADBE0002`) — Camera orbit limits (elevation min/max, radius min) for restricting the view. Implemented in `extensions/cc/safe-orbit-camera-adobe.h` and `safe-orbit-camera-adobe.cc`. See the main [README](../README.md) for attributes and defaults.

- **SPZ_ADOBE_coordinate_system** (`0xADBE0003`) — Records the coordinate system in which the Gaussian data is physically stored in the file. Implemented in `extensions/cc/coordinate-system-adobe.h` and `coordinate-system-adobe.cc`.

  **Payload:** one `uint32_t` — the `CoordinateSystem` enum value (4 bytes). Valid range: 1–16 (see `CoordinateSystem` enum). Value 0 (`UNSPECIFIED`) is treated as absent and will produce a warning; it is not a valid way to express "use the default".

  **This extension is a descriptor, not an instruction.** It labels what coordinate system the packed data is already in. Readers must not apply an additional conversion based on it — `unpackGaussians` uses it internally and performs the full conversion to `UnpackOptions::to` automatically. Manually applying a conversion on top will double-transform the data.

  **For writers:** attach this extension before saving to store data in a coordinate system other than the default RUB. `packGaussians` will then convert `PackOptions::from → extension.coordinateSystem` instead of `PackOptions::from → RUB`. The library never creates this extension automatically.

  **For readers:** if you use `unpackGaussians` / `loadSpz`, no action is needed — the conversion is handled for you. Only inspect this extension directly if you are processing raw packed data outside of those functions.

  **Compatibility note:** a non-extension build ignores this extension and always assumes RUB, logging a warning at both load and unpack time. If the file was stored in a different coordinate system, the resulting data will be incorrect. See "Builds without extension support" above.

  ```python
  ext = spz.SpzExtensionCoordinateSystemAdobe()
  ext.coordinate_system = spz.RDF
  cloud.extensions = [ext]
  # pack with from_coord=RUB → stored as RDF; load with to_coord=RUB → converted back automatically
  ```

- **SPZ_NIANTIC_georeference** (`0x4E530001`) — Records a similarity transform that places the asset in a geocentric CRS. Implemented in `extensions/cc/georeference-niantic.h` and `georeference-niantic.cc`.

  **Payload** (little-endian, byte-packed, 76 bytes fixed + variable `crs`):

  | Offset | Size | Field | Description |
  |--------|------|-------|-------------|
  | 0 | 1 | `ext_version` | `uint8_t`, must be 1 |
  | 1 | 1 | `flags` | `uint8_t`, bit 0 `has_epoch`, bits 1-2 `crs_encoding`, bits 3-7 reserved |
  | 2 | 2 | `crs_length` | `uint16_t`, bytes of the `crs` string |
  | 4 | 24 | `origin[3]` | `double`, meters |
  | 28 | 32 | `rotation[4]` | `double`, unit quaternion x, y, z, w |
  | 60 | 8 | `scale` | `double`, uniform scale |
  | 68 | 8 | `epoch` | `double`, coordinate epoch as a decimal year; writers write NaN when `has_epoch` = 0 |
  | 76 | `crs_length` | `crs` | Target CRS, **not** null-terminated; syntax given by `crs_encoding` |

  **Transform:** `p_crs = scale * (q * p_local * q⁻¹) + origin`, where `p_local` is a position in RUB, whatever frame the data is stored in. `q` is `rotation`, a unit Hamilton quaternion (x, y, z, w) applied as an active rotation from RUB axes to the target CRS axes (e.g. ECEF, not a local east-north-up frame). `scale` is uniform and positive; `origin` is in meters. The library never applies the transform to the Gaussian data.

  **Target CRS:** `crs` is the only CRS in the record. The target must be a geocentric CRS: a geodetic CRS with X, Y, Z axes in meters (ECEF on Earth). Heights derived from `p_crs` are ellipsoidal. `crs_encoding` gives the syntax of `crs`:

  - `0` = an `AUTHORITY:CODE` naming a CRS in the PROJ database, e.g. `EPSG:4978` (WGS 84 ECEF) or `EPSG:9988` (ITRF2020). Write it exactly as the database lists it (`EPSG:4978`, not `epsg:4978`), so readers can compare codes as strings.
  - `1` = an inline PROJJSON `GeodeticCRS` with a Cartesian coordinate system, for a target without a code — e.g. the Moon, using the datum of `IAU_2015:30100`.

  Writers use a code when the target has one. A reader that cannot resolve `crs` may still apply the transform, but must not assume any body or datum. The library performs no registry lookup, so conformance of the CRS type is the writer's responsibility.

  **Epoch:** `epoch` is the coordinate epoch of `p_crs` (OGC 18-010r11 §16): the decimal year at which the coordinates are valid in a dynamic target frame (e.g. ITRF2020 at 2026.7), not the frame's reference epoch. `crs` carries no coordinate epoch: no `@epoch` code suffix and no PROJJSON `CoordinateMetadata`. Writers set `epoch` for dynamic targets and may set it for `EPSG:4978`, whose WGS 84 datum ensemble has only dynamic members; that epoch does not identify a WGS 84 realization, so use a realization-specific CRS when that matters. Readers ignore it for static targets.

  **Validation on load** (a rejected payload is skipped with a warning and the core data still loads). The payload is rejected unless:

  - `ext_version` is 1;
  - `crs_encoding` is `0` or `1`;
  - `crs_length` does not exceed the remaining payload bytes;
  - `crs` is non-empty and, for `crs_encoding` = 0, printable ASCII (no spaces, no `@`) with exactly one `:` and both parts non-empty, or, for `crs_encoding` = 1, well-formed UTF-8 (it is not parsed). Unknown authorities and codes are accepted as written;
  - every `origin` component is finite;
  - `rotation` satisfies `abs(sqrt(x² + y² + z² + w²) − 1) <= 1e-6`;
  - `scale` is finite and greater than zero;
  - `epoch` is finite when `has_epoch` is set.

  Reserved `flags` bits are ignored on read; writers write them as 0.

  **Validation on write:** a `crs` that fails the rules above or exceeds 65,535 bytes is omitted with a warning, so readers reject the record.

  **For writers:** attach this extension before saving, with the transform defined from RUB positions. The library never creates this extension automatically.

  **For readers:** loaded positions are always local; the library never applies the transform. Load with `to_coord` = RUB and apply the transform yourself in double precision.

  **Compatibility note:** a non-extension build skips this extension with a warning. The Gaussian data is unaffected; only the georeference is lost.

  ```python
  ext = spz.SpzExtensionGeoreferenceNiantic()
  ext.crs = "EPSG:4978"                       # WGS 84 ECEF; PROJJSON + crs_encoding for other bodies
  ext.origin = [4194304.5, -555555.25, 4713930.125]
  ext.rotation = [0.0, 0.0, 0.0, 1.0]         # x, y, z, w
  ext.scale = 1.0
  cloud.extensions = [ext]

  # Reading: load in RUB, then apply the transform to the local positions
  from scipy.spatial.transform import Rotation  # from_quat expects x, y, z, w
  opts = spz.UnpackOptions()
  opts.to_coord = spz.CoordinateSystem.RUB
  loaded = spz.load_spz("georeferenced.spz", opts)
  geo = next(e for e in loaded.extensions
             if e.extension_type == spz.SpzExtensionType.SPZ_NIANTIC_georeference)
  p_local = loaded.positions.reshape(-1, 3).astype(np.float64)
  p_ecef = geo.scale * Rotation.from_quat(geo.rotation).apply(p_local) + np.array(geo.origin)
  ```
