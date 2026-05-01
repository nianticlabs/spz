# spz

`.spz` is a file format for compressed 3D gaussian splats. This directory contains a C++ library
for saving and loading data in the .spz format.

spz encoded splats are typically around 10x smaller than the corresponding .ply files,
with minimal visual differences between the two.

## Internals

### Coordinate System

SPZ stores data internally in an RUB coordinate system following the OpenGL and three.js
convention. This differs from other data formats such as PLY (which typically uses RDF), GLB (which
typically uses LUF), or Unity (which typically uses RUF). To aid with coordinate system conversions,
callers should specify the coordinate system their Gaussian Cloud data is represented in when saving
and what coordinate system their rendering system uses when loading. These are specified in the
PackOptions and UnpackOptions respectively.  If the coordinate system is `UNSPECIFIED`, data will
be saved and loaded without conversion, which may harm interoperability.

There are 16 named coordinate systems organised in two families:

- **Standard family** — letters denote the X, Y, Z axis directions (e.g. `RUB` = Right, Up, Back):
  `LDB`, `RDB`, `LUB`, `RUB`, `LDF`, `RDF`, `LUF`, `RUF`
- **Rotated family** — axes are permuted by a 90-degree rotation about X (Y and Z are swapped):
  `LBD`, `RBD`, `LBU`, `RBU`, `LFD`, `RFD`, `LFU`, `RFU`

Converting between families applies a 90-degree rotation about the X axis in addition to any axis
flips, which also rotates spherical-harmonics coefficients via the appropriate Wigner D-matrix.

## Implementations

### C++

Requires `libz` and `libzstd` as dependent libraries, otherwise the code is completely self-contained.
A CMake build system is provided for convenience.

**Note:** If `libz` or `libzstd` are not found on the system, the CMake build system will automatically download them using FetchContent. This ensures consistent builds across different environments.

### Typescript

To build the Typescript interface through Web-Assembly (WASM), an Emscripten environment needs to be setup before compilation. One may install the Emscripten SDK following the instructions [here](https://emscripten.org/docs/getting_started/downloads.html).

Under an emscripten environment, the Makefile can be generated through:
```
emcmake cmake -B build-wasm .
```
Then one can build through
```
cmake --build build-wasm
```
The package will be built and installed into the `dist` folder.


## API

### C++

```
bool saveSpz(
   const GaussianCloud &gaussians, const PackOptions &options, std::vector<uint8_t> *output);
```

Converts a cloud of Gaussians to `.spz` format and writes the result to a vector of bytes.

   - `gaussians`: The Gaussians to save
   - `options`: Flags that control the packing behavior, including SH quantization parameters.
   - `output`: A vector that will be populated with bytes encoded in .spz format
   - Returns true on success and false on failure.

---

```
bool saveSpz(
   const GaussianCloud &gaussians, const PackOptions &options, const std::string
&filename);
```

Saves a cloud of Gaussians in `.spz` format to a file

   - `gaussians`: The Gaussians to save
   - `options`: Flags that control the packing behavior, including SH quantization parameters.
   - `filename`: The path to the file to save to.
   - Returns true on success and false on failure.

---

```
GaussianCloud loadSpz(const std::vector<uint8_t> &data, const UnpackOptions &opti
ons);
```

Loads a cloud of Gaussians from bytes in `.spz` format.

   - `data`: A vector containing the encoded spz data
   - `options`: Flags that control the unpacking behavior.
   - Returns a `GaussianCloud` decoded from the vector. In case of an error, this will return
     a result with no gaussians

---

```
GaussianCloud loadSpz(const std::string &filename, const UnpackOptions &options);
```

Loads a cloud of Gaussians from a file in `.spz` format.

   - `filename`: The path to the file to load from.
   - `options`: Flags that control the unpacking behavior.
   - Returns a `GaussianCloud` decoded from the file. In case of an error, this will return
     a result with no gaussians

### Typescript

Check [src/emscripten/spz.d.ts.in](src/emscripten/spz.d.ts.in) (source) or `dist/spz.d.ts` (after building with CMake) for the TypeScript interface. Since the Emscripten and Javascript memory are separately handled, we only expose limited functionalities for the Typescript interface.

### PackOptions

The `PackOptions` struct supports the following fields:

- `from`: Source coordinate system (default: `UNSPECIFIED`)
- `version`: Version of the packed format (default: `4`)
- `sh1Bits`: Number of quantization bits for SH degree 1 coefficients (default: 5, range: 1-8)
- `shRestBits`: Number of quantization bits for SH degree 2+ coefficients (default: 4, range: 1-8)

## File Format

Version 4 uses ZSTD compression with a 32-byte plaintext header. Each attribute stream is compressed independently and in parallel, improving both compression and decompression throughput. Files from versions 1–3 (gzip-compressed single-stream format) can still be read. There is no maximum point count enforced during loading.

### File layout (version 4)

```
Bytes  0–31:              NgspFileHeader (32 bytes, plaintext)
Bytes  32..(tbo-1):       Extension ILV records, if flags & 0x2  (variable)
Bytes  tbo..(tbo+N*16-1): TOC — N × [compressedSize u64, uncompressedSize u64]
Bytes  (tbo+N*16)..end:   N independent ZSTD-compressed attribute streams
```

`tbo` = `tocByteOffset` from the header. When no extensions are present, `tbo = 32`.

### Header (version 4)

```c
struct NgspFileHeader {        // 32 bytes, all fields little-endian
  uint32_t magic;              // 0x5053474e ("NGSP")
  uint32_t version;            // 4
  uint32_t numPoints;
  uint8_t  shDegree;
  uint8_t  fractionalBits;
  uint8_t  flags;
  uint8_t  numStreams;         // number of ZSTD-compressed attribute streams (typically 6)
  uint32_t tocByteOffset;     // byte offset from file start to the TOC
  uint8_t  reserved[12];      // zero, reserved for future use
};
```

All values are little-endian.

1. **magic**: Always `0x5053474e` (bytes `N`, `G`, `S`, `P` in file order).
2. **version**: 4 (current). Versions 1–3 use the legacy gzip format (read-only).
3. **numPoints**: The number of gaussians.
4. **shDegree**: The degree of spherical harmonics. Must be between 0 and 4 (inclusive).
5. **fractionalBits**: The number of bits used for the fractional part of fixed-point coordinates.
6. **flags**: Bit field.
   - `0x1`: splat was trained with [antialiasing](https://niujinshuchong.github.io/mip-splatting/).
   - `0x2`: the header zone contains vendor-specific extension records (see [Extensions](#extensions)).
7. **numStreams**: Number of ZSTD-compressed attribute streams following the TOC. Typically 6 (positions, alphas, colors, scales, rotations, spherical harmonics).
8. **tocByteOffset**: Byte offset from the start of the file to the Table of Contents. Everything before `tocByteOffset` is plaintext.
9. **reserved**: Must be zero.

### Format detection (load path)

```
bytes[0..3] == "NGSP"  →  version 4 ZSTD format
bytes[0..1] == 0x1f 0x8b  →  legacy gzip format (versions 1–3, read-only)
```

### Positions

Positions are represented as `(x, y, z)` coordinates, each as a 24-bit fixed point signed integer.
The number of fractional bits is determined by the `fractionalBits` field in the header.

### Scales

Scales are represented as `(x, y, z)` components, each represented as an 8-bit log-encoded integer.

### Rotation

In version 3 and 4, rotations are represented as the smallest three components of the normalized rotation quaternion, for optimal rotation accuracy.
The largest component can be derived from the others and is not stored. Its index is stored on 2 bits
and each of the smallest three components is encoded as a 10-bit signed integer.

In version 2, rotations are represented as the `(x, y, z)` components of the normalized rotation quaternion. The
`w` component can be derived from the others and is not stored. Each component is encoded as an
8-bit signed integer.

### Alphas

Alphas are represented as 8-bit unsigned integers.

### Colors

Colors are stored as `(r, g, b)` values, where each color component is represented as an
unsigned 8-bit integer.

### Spherical Harmonics

Depending on the degree of spherical harmonics for the splat, this can contain 0 (for degree 0),
9 (for degree 1), 24 (for degree 2), 45 (for degree 3), or 72 (for degree 4) coefficients per gaussian.

The coefficients for a gaussian are organized such that the color channel is the inner (faster
varying) axis, and the coefficient is the outer (slower varying) axis, i.e. for degree 1,
the order of the 9 values is:

```
sh1n1_r, sh1n1_g, sh1n1_b, sh10_r, sh10_g, sh10_b, sh1p1_r, sh1p1_g, sh1p1_b
```

Each coefficient is represented as an 8-bit signed integer.

**Quantization:**

SPZ supports configurable spherical harmonics quantization. By default, a fixed-precision quantization is adopted:
- Assumes SH coefficients are in the [-1, 1] range
- Fixed 5 bits of precision for degree 1 and 4 bits for degrees 2, 3, and 4
- Maintained for backward compatibility

The quantization precision can be configured via `PackOptions`:
- `sh1Bits`: Number of quantization bits for SH degree 1 coefficients (default: 5, range: 1-8)
- `shRestBits`: Number of quantization bits for SH degree 2+ coefficients (default: 4, range: 1-8)

**Note:** Quantization bits are only used during packing to reduce information entropy for better compression. The unpacking process does not need to know the exact quantization bits, as decompression fills zero bits for quantized data automatically.

This allows users to trade off between file size and quality. The library maintains full backward compatibility with default quantization settings.

### Extensions

SPZ supports vendor-specific extensions (e.g. camera limits) so multiple vendors can coexist in the same file. In version 4, extensions are stored as ILV (type + length + value) records in the plaintext header zone between the fixed header and the TOC — parsers that don't recognize an extension type skip it by length. For the extension format and how to add or use extensions, see [extensions/README.md](extensions/README.md).

### Camera Orbit Limitation

With extension `SPZ_ADOBE_safe_orbit_camera`, SPZ supports storing camera limits which can be used to restrict the view in a render. This extension includes:

**Attributes**

This extension has the following attributes and default values:

```
  float safeOrbitElevationMin = 0.0f;  // Minimum elevation for safe orbit (radians)
  float safeOrbitElevationMax = 0.0f;  // Maximum elevation for safe orbit (radians)
  float safeOrbitRadiusMin = 0.0f;     // Minimum radius for safe orbit
```


## Python Bindings

The SPZ library provides Python bindings built with [nanobind](https://nanobind.readthedocs.io/) that offer a convenient interface for loading, manipulating, and saving 3D Gaussian splats from Python.

### Installation
```bash
git clone https://github.com/nianticlabs/spz.git
cd spz
pip install .
```

Please see src/python/README.md for more details and usage examples

## Web Utility

A browser-based utility is available at [nianticlabs.github.io/spz](https://nianticlabs.github.io/spz). It uses the latest WASM build and allows you to inspect file metadata for `.spz` and `.ply` files and convert between the two formats without installing anything. As it runs entirely in the browser, it is subject to browser memory limits and is intended for files up to around 2 GB. For larger files, clone the repository and use the command-line tools directly.
