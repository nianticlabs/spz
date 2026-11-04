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
   Extensions are part of `PackedGaussians` and `GaussianCloud`. When you call `loadSpz` / `loadSpzPacked` or `saveSpz`, the extension list is read/written automatically. For version 4 and above, `saveSpz` writes extension records into the plaintext header zone, and any bulk data they carry into extension frames (see [Storing bulk data in extension frames](#storing-bulk-data-in-extension-frames)).

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

When loading PLY files, the loader walks the non-vertex elements in file order. Each extension that supports PLY declares the element names it uses in a static member (e.g. `SpzExtensionSafeOrbitCameraAdobe::kRequiredPlyElementNames`) and is registered with an exemplar in `getPlyExtensionRegistry()` in `extensions/cc/splat-extensions.cc`. When an element belongs to a registered extension, `readExtensionsFromPly` hands that extension the run of consecutive elements it declares; any other element is skipped by its size (element count times the summed property sizes: `char`/`uchar` 1, `short`/`ushort` 2, `int`/`uint`/`float` 4, `double` 8). After each run the reader moves to the end of it, so an extension that rejects its data never shifts what follows, and unknown elements (e.g. `up_axis`, `chirality`) may sit before, between or after known ones. The union of all registered element names answers **`isKnownPlyExtensionElement(const std::string& elementName)`**. To add PLY support for a new extension, implement the three abstract PLY methods on your extension struct and register it (see step 6 under "Adding a new extension").

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
   - Optional: `bool isValid() const override` — return false, logging the reason, if `write()` would emit a payload your `read()` rejects; `saveSpz` then fails. Defaults to true.
   - `static std::optional<SpzExtensionBasePtr> read(std::istream& is)` — read payload from `is` (caller provides a stream of exactly your payload bytes).
   - `static SpzExtensionType type()` — return your enum value.
   - The three PLY methods (required because they are pure virtual on the base): `tryReadFromPly(std::istream& in, const std::vector<PlyExtraElement>& elements) const override`, `writePlyHeader(std::ostream& out) const override`, `writePlyData(std::ostream& out) const override`. If your extension does not support PLY, implement them to return `std::nullopt` / no-op as appropriate; if it does, see step 6.

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
   - **Implement the three PLY methods** — `tryReadFromPly(in, elements)`: `elements` is the run of consecutive elements (name, count, bytes per element) from your `kRequiredPlyElementNames`, in file order, and `in` is at the start of their data; check the layout, read, and return a new instance (or `std::nullopt` on failure; the caller then skips the rest of the run). `writePlyHeader(out)`: write the PLY `element` / `property` lines for your data. `writePlyData(out)`: write the binary data for your extension, in the same order. `saveSplatToPly` validates extensions first and converts them from `PackOptions::from` to RDF like the vertex data.
   - **Register in the PLY registry** — In `extensions/cc/splat-extensions.cc`, inside `getPlyExtensionRegistry()`, add an entry: `{&MyExtension::kRequiredPlyElementNames, std::make_shared<MyExtension>()}`. You do **not** edit `readExtensionsFromPly`, `writeExtensionsToPlyHeader`, or `writeExtensionsToPlyData`; they dispatch to your virtual methods.

7. **Bindings**  
   Expose the new type and struct in Python (`extensions/python/splat-extensions.cc`) and Emscripten (`extensions/emscripten/splat-extensions.cc`, plus `.d.ts.in` if needed) so callers can construct and read your extension from the cloud’s `extensions` list. Each binding file must `#include` your extension’s header (e.g. `#include "safe-orbit-camera-adobe.h"`) to get the full type definition.

### Storing bulk data in extension frames

Records live in the plaintext zone before the TOC, which is capped at 4 GiB and must be read before any stream. Extensions that carry per-point or other bulk data should store it in extension frames instead: zstd skippable frames, each wrapping one zstd frame, appended to the last core stream (see [Extension frames](../README.md#extension-frames)). The record keeps only a small descriptor and the list of its extension frames.

Optional capabilities are separate interfaces, declared in `splat-extensions.h`, that an extension inherits alongside `SpzExtensionBase`. The library detects them with `dynamic_cast`, so `SpzExtensionBase`, and every extension built against it, is unchanged.

| Interface | What the extension provides | What the library does |
| --- | --- | --- |
| `SpzFramedExtension` | `planFrameWrites(numPoints)`: an `SpzFrameWritePlan` holding the bulk data and a `writeRecord(os, frames)` that writes the record once its extension frames are placed (`std::nullopt` fails the save) | Splits the data into extension frames of at most 16 MiB uncompressed, compresses them in parallel with the core streams, and passes one `SpzFrameLocation` per extension frame to `writeRecord`. `writeAllExtensions` skips framed records, since they can't be written without extension frames. |
| `SpzFrameReadPlan` (a struct returned by the record's parser) | `type`, the extension `frames`, an optional uniqueness `key`, and two functions: `dataBytes(numPoints)` and `resolve(numPoints, data)` | Drops records sharing a key, validates the extension frames against the file and the decoded size before allocating, decompresses them in parallel with the core streams, and calls `resolve` with the concatenated data. A record that fails any step is dropped; the core data still loads. |
| `SpzKeyedExtension` | `uniquenessKey()` | `saveSpz` fails if two records of the same type share a key; loading drops every record that does. |
| `SpzCoordinateDependentExtension` | `withConvertedCoordinates(c)`: a converted copy, or `nullptr` if unaffected | Called whenever the cloud is converted between coordinate systems. Records are shared between clouds, so never modify `*this`. |

Extension frames exist only in version 4 files: `saveSpz` refuses older versions when a record has frame data.

**Why extension frames are listed by their record, not the TOC.** Released readers require `numStreams` to equal the number of core streams, and the TOC sizes to add up to the end of the file, so an extra TOC entry would make them reject the file. A record's frame table is invisible to them, because unknown records are skipped by length, and the extension frames still count toward the file length through the last core stream's `compressedSize`. zstd would skip them in any stream; the last one is used because the core streams then stay contiguous at the front of the file, with all extension data at its end. If a future v5 bumps the format version anyway, moving extension frames into real TOC entries would be the cleaner end state.

## Built-in extensions

- **SPZ_NIANTIC_custom_attribute** (`0x4E530002`) — A named per-point attribute (normals, material parameters, embeddings, ...); a cloud may carry any number of them. Implemented in `extensions/cc/custom-attribute-niantic.h` and `custom-attribute-niantic.cc`. The values are stored in [extension frames](#storing-bulk-data-in-extension-frames), so attribute-bearing files load in every SPZ v4 reader; readers without this extension get the core cloud and skip the attribute.

  `dim` is the number of floats per point (1–65535, so embeddings of any width fit). `semantic` says whether coordinate conversions apply:

  | `semantic` | Value | `dim` | Coordinate conversion | Typical use |
  | --- | --- | --- | --- | --- |
  | `CUSTOM` | 0 | Any | Never converted (like glTF's `_`-prefixed attributes) | Material parameters, embeddings, IDs |
  | `SPATIAL` | 1 | Multiple of 3 | Each xyz triple rotated and flipped like positions | Normals, directions, velocities |

  The `name` is a free-form label used to look the attribute up (`"roughness"`, a PLY property name, …); it has no effect on behavior. Names are unique within a cloud: `saveSpz` fails on a duplicate, and loading drops every attribute that shares a name. Readers drop attributes whose `semantic` they don't know, since they can't know how those convert.

  **Payload:**

  ```
  u16 nameLength, u8 name[nameLength]   // printable ASCII (0x20-0x7E), 1-65535 bytes
  u16 dim                               // floats per point, >= 1
  u8  semantic                          // 0 CUSTOM, 1 SPATIAL (dim % 3 == 0)
  u8  bits                              // 32 = float32; 1-16 = quantized
  f32 min[dim], f32 max[dim]            // only when quantized
  u32 frameCount
  frameCount x { u64 frameOffset; u64 compressedSize; u64 uncompressedSize; }   // the frame table: one
                                                                                 // entry per extension frame
  ```

  Joining the extension frames' decompressed payloads in order gives the values component-major: every point's first component, then every point's second, and so on. float32 components are 4 bytes, little-endian. Quantized components store level `round((v - min) / (max - min) * (2^bits - 1))`, rounding halves away from zero, or 0 when `min == max`; decoded in double precision as `min + level * ((max - min) / (2^bits - 1))`. A level takes `ceil(bits / 8)` bytes, stored as byte planes: a component's data is every point's least significant level byte, then (for `bits` 9–16) every point's most significant byte. Splitting the bytes this way compresses better than storing them together; with `bits` ≤ 8 there is a single plane, one byte per point.

  Readers MUST drop the attribute (the core cloud still loads) when:
  - `nameLength` is 0, or the name has a byte outside printable ASCII (0x20-0x7E);
  - `semantic` is unknown, or `SPATIAL` with `dim` not a multiple of 3, or `dim` is 0;
  - `bits` is not 32 or 1–16, a quantized range isn't finite with `min <= max`, or a decoded level exceeds `2^bits - 1`;
  - the extension frames' decompressed sizes don't add up to `numPoints * dim * (bits == 32 ? 4 : ceil(bits / 8))`;
  - another custom attribute in the file with a well-formed descriptor (valid name, `dim`, `semantic` and `bits`) has the same name, compared byte for byte without Unicode normalization. This applies even when that attribute's extension frames or values are invalid; an attribute whose descriptor is malformed is dropped on its own and doesn't count.

  Readers MUST ignore bytes after the frame table: later versions may append fields there.

  An attribute of any size is supported; large ones simply have more extension frames. Readers that can't hold the decoded attribute in memory (e.g. over 4 GiB on wasm32) drop it and still load the core cloud.

  `SPATIAL` values are stored in the same coordinate system as the file's positions (RUB, or the system named by `SPZ_ADOBE_coordinate_system`) and are converted with them on load. The conversion treats each triple as an ordinary (polar) vector, which suits normals, directions and velocities; it would give pseudovectors such as angular velocity the wrong sign under mirroring conversions.

  Readers that predate this extension, and builds without extension support, load the core cloud and log a warning that an unknown extension was skipped and "data may be incorrect". For custom attributes that warning is harmless: they never change how the core data is packed. PLY can't carry custom attributes yet; saving to PLY drops them with a warning.

  ```python
  normal = spz.SpzExtensionCustomAttributeNiantic()
  normal.name = "normal"
  normal.dim = 3
  normal.semantic = spz.SpzAttributeSemantic.SPATIAL  # converted with coordinates
  normal.bits = 32
  normal.values = normals                             # (num_points, 3) floats

  embedding = spz.SpzExtensionCustomAttributeNiantic()
  embedding.name = "clip_embedding"                   # semantic defaults to CUSTOM: never converted
  embedding.dim = 512
  embedding.bits = 8
  embedding.values = clip_features                    # num_points * 512 floats

  cloud.extensions = [normal, embedding]

  options = spz.UnpackOptions()
  options.custom_attributes = ["normal"]              # attributes load only when asked for
  loaded = spz.load_spz("scene.spz", options)
  normals = spz.find_custom_attribute(loaded, "normal").values.reshape(-1, 3)
  ```

  In C++, `spz::findCustomAttribute(cloud.extensions, "normal")` does the same; in JavaScript, attributes expose `getValues()` and `setValues()`.

  Custom attributes are opt-in on load: `UnpackOptions::customAttributes` (`custom_attributes` in Python, `customAttributes` in JavaScript) names the ones to decode, and its default, an empty list, decodes none, just like a build without extension support. Attributes cost memory in proportion to the point count, and most loaders (viewers) only render the core cloud. Unlisted attributes are skipped before their extension frames are decompressed, so they cost nothing. `std::nullopt` / `None` / `null` decodes them all; **tools that load and re-save a file must set it, or the attributes are lost**.

  ```js
  const viewer = spz.loadSpzFromBuffer(bytes, { to: spz.CoordinateSystem.RUB });  // no attributes
  const all = spz.loadSpzFromBuffer(bytes, { to: spz.CoordinateSystem.RUB, customAttributes: null });
  ```


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

  **Epoch:** `epoch` is the coordinate epoch of `p_crs` (OGC 18-010r11 §16): the decimal year at which the coordinates are valid in a dynamic target frame (e.g. ITRF2020 at 2026.7), not the frame's reference epoch. `crs` carries no coordinate epoch: no `@epoch` code suffix and no PROJJSON `CoordinateMetadata`. Writers set `epoch` for dynamic targets and should not set it for `EPSG:4978`: an epoch does not identify a realization of the WGS 84 datum ensemble, so use a realization-specific CRS (e.g. `EPSG:9988`) when the epoch matters. Readers ignore it for static targets.

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

  **Validation on save:** a `rotation` within 1e-3 of unit norm is normalized. If the record then breaks any rule above, including a `crs` over 65,535 bytes, `saveSpz` logs the reason and returns false without saving.

  **For writers:** attach this extension before saving, with the transform defined from RUB positions. The library never creates this extension automatically.

  **For readers:** loaded positions are always local; the library never applies the transform. Load with `to_coord` = RUB and apply the transform yourself in double precision.

  **Compatibility note:** a non-extension build skips this extension with a warning. The Gaussian data is unaffected; only the georeference is lost.

  ```python
  import numpy as np

  ext = spz.SpzExtensionGeoreferenceNiantic()
  ext.crs = "EPSG:9988"                       # ITRF2020 ECEF; EPSG:4978 for WGS 84; PROJJSON + crs_encoding for other bodies
  ext.epoch = 2026.7                          # coordinate epoch, decimal year
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

- **SPZ_ADOBE_environment_map** (`0xADBE0005`) — An equirectangular RGB environment map with the frame it is looked up in, e.g. a baked reflection environment. At most one per cloud. Implemented in `extensions/cc/environment-map-adobe.h` and `environment-map-adobe.cc`.

  | Field | Meaning |
  |---|---|
  | `width`, `height` | Map size in pixels, both at least 1 |
  | `origin` | Lookup origin; converted like positions by `PackOptions::from`, `UnpackOptions::to` and `convertCoordinates` |
  | `up`, `e1`, `e2` | Unit directions, converted the same way. Pixel column = azimuth about `up` measured from `e1` toward `e2`, -pi..pi → 0..width; pixel row = angle from `up`, 0..pi → 0..height. Because the frame follows the splats, lookups stay correct in any coordinate system |
  | `metadata` | Free-form UTF-8 (JSON by convention, e.g. renderer settings), never converted |
  | `pixels` | `height * width * 3` floats, row-major, RGB interleaved, linear-in-file values in the model's color space |

  **Payload** (little-endian):

  ```
  u32 width, u32 height
  f32 origin[3], f32 up[3], f32 e1[3], f32 e2[3]
  u32 metadataLength, u8 metadata[metadataLength]
  zstd frame, to the end of the payload: height * width * 3 float16 (IEEE half bits), row-major, RGB interleaved
  ```

  Pixels are float32 in memory and float16 in the SPZ record, so float16 inputs round-trip exactly. Saving fails on a zero width or height, a pixel count other than `width * height * 3`, non-finite pixels or vectors, pixels beyond ±65504, or a second environment map. On load, a malformed record is dropped with an error and the core splat data still loads; a file with more than one environment map loses all of them, like any extension whose records share a uniqueness key.

  **PLY form** (after `vertex`, these three elements in this order; binary little-endian; pixels stay float32):

  ```
  element environment 1
  property uint environment_width
  property uint environment_height
  property float environment_origin_x      // then _y, _z; likewise environment_up_*, environment_e1_*, environment_e2_*
  ...
  element environment_map <width*height>
  property float environment_map_r
  property float environment_map_g
  property float environment_map_b
  element environment_metadata <metadataLength>   // count 0 allowed
  property uchar environment_metadata
  ```

  ```python
  env = np.load("reflection_env.npy")  # (H, W, 3) float16
  cloud.extensions = [spz.SpzExtensionEnvironmentMapAdobe(
      env.shape[1], env.shape[0], env.astype(np.float32).ravel(),
      origin=(0, 0, 0), up=(0, 1, 0), e1=(1, 0, 0), e2=(0, 0, 1), metadata='{"mip_levels": 6}')]
  spz.save_spz(cloud, spz.PackOptions(), "out.spz")  # or spz.save_splat_to_ply(...)

  loaded = spz.load_spz("out.spz")
  (env_map,) = [e for e in loaded.extensions if isinstance(e, spz.SpzExtensionEnvironmentMapAdobe)]
  image = env_map.pixels.reshape(env_map.height, env_map.width, 3)
  ```
