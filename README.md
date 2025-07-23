# spz

`.spz` is a file format for compressed 3D gaussian splats. This directory contains a C++ library
for saving and loading data in the .spz format.

spz encoded splats are typically around 10x smaller than the corresponding .ply files,
with minimal visual differences between the two.

### Discord
Join the [Scaniverse Discord](https://discord.gg/xKtCxvEmxa) to discuss SPZ-related topics with Niantic developers and community contributors. 

## Internals

### Coordinate System

SPZ stores data internally in an RUB coordinate system following the OpenGL and three.js
convention. This differs from other data formats such as PLY (which typically uses RDF), GLB (which
typically uses LUF), or Unity (which typically uses RUF). To aid with coordinate system conversions,
callers should specify the coordinate system their Gaussian Cloud data is represented in when saving
and what coordinate system their rendering system uses when loading. These are specified in the
PackOptions and UnpackOptions respectively.  If the coordinate system is `UNSPECIFIED`, data will
be saved and loaded without conversion, which may harm interoperability.

### Spherical Harmonics Quantization

Starting with version 3, SPZ supports configurable spherical harmonics quantization with the following improvements:

1. **Adaptive Range Scaling**: Instead of assuming SH coefficients are in the [-1, 1] range, the format now computes the actual min/max values from the data and uses the full quantization range for optimal precision.

2. **Configurable Bit Precision**: The number of quantization bits for SH degree 1 and higher degree coefficients can be configured (1-8 bits), allowing users to trade off between file size and quality.

3. **Backward Compatibility**: The library maintains full backward compatibility with versions 1 and 2, automatically detecting and using legacy quantization methods for older files.

## Implementations

### C++

Requires `libz` as the only dependent library, otherwise the code is completely self-contained.
A CMake build system is provided for convenience.

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

To publish the package to the NPM repository, one may login in with
```
npm login --registry https://artifactory.corp.adobe.com/artifactory/api/npm/npm-adobe-release-local
```
Then publish through
```
npm run publishPackage
```


## API

### C++

```
bool saveSpz(
   const GaussianCloud &gaussians, PackOptions &options, std::vector<uint8_t> *output);
```

Converts a cloud of Gaussians in `.spz` format to a vector of bytes.

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

Check [src/emscripten/spz.d.ts](src/emscripten/spz.d.ts) for the Typescript interface. Since the Emscripten and Javascript memory are separately handled, we only expose limited functionalities for the Typescript interface.

### PackOptions

The `PackOptions` struct supports the following fields:

- `from`: Source coordinate system (default: `UNSPECIFIED`)
- `sh1Bits`: Number of quantization bits for SH degree 1 coefficients (default: 5, range: 1-8)
- `shRestBits`: Number of quantization bits for SH degree 2+ coefficients (default: 4, range: 1-8)
- `hasSafeOrbit`: Whether safe orbit camera data is present (default: false)
- `safeOrbitElevationMin`: Minimum elevation for safe orbit in radians (default: 0.0)
- `safeOrbitElevationMax`: Maximum elevation for safe orbit in radians (default: 0.0)
- `safeOrbitRadiusMin`: Minimum radius for safe orbit (default: 0.0)

## File Format

The .spz format is a gzipped stream of data consisting of a variable-size header followed by the
gaussian data. This data is organized by attribute in the following order: positions,
alphas, colors, scales, rotations, spherical harmonics.

### Header

**Version 4 (current):**
```c
struct PackedGaussiansHeader {
  uint32_t magic;
  uint32_t version;
  uint32_t numPoints;
  uint8_t shDegree;
  uint8_t fractionalBits;
  uint8_t flags;
  uint8_t reserved;
  uint8_t sh1Bits;                   // Bits for SH degree 1 coefficients
  uint8_t shRestBits;                // Bits for SH degree 2+ coefficients  
  float shMin;                       // Minimum SH coefficient value for quantization
  float shMax;                       // Maximum SH coefficient value for quantization
  uint8_t hasSafeOrbit;              // Whether safe orbit data is present
  uint8_t padding[1];                // Padding for alignment
  float safeOrbitElevationMin;       // Minimum elevation for safe orbit (radians)
  float safeOrbitElevationMax;       // Maximum elevation for safe orbit (radians)
  float safeOrbitRadiusMin;          // Minimum radius for safe orbit
};
```

**Version 3 (legacy):**
```c
struct PackedGaussiansHeader {
  uint32_t magic;
  uint32_t version;
  uint32_t numPoints;
  uint8_t shDegree;
  uint8_t fractionalBits;
  uint8_t flags;
  uint8_t sh1Bits;        // Bits for SH degree 1 coefficients
  uint8_t shRestBits;     // Bits for SH degree 2+ coefficients  
  uint8_t reserved;
  float shMin;            // Minimum SH coefficient value for quantization
  float shMax;            // Maximum SH coefficient value for quantization
};
```

**Version 2 (legacy):**
```c
struct PackedGaussiansHeader {
  uint32_t magic;
  uint32_t version;
  uint32_t numPoints;
  uint8_t shDegree;
  uint8_t fractionalBits;
  uint8_t flags;
  uint8_t reserved;
};
```

All values are little-endian.

1. **magic**: This is always 0x5053474e
2. **version**: Currently supported versions are 1, 2, 3, and 4
3. **numPoints**: The number of gaussians
4. **shDegree**: The degree of spherical harmonics. This must be between 0 and 4 (inclusive).
5. **fractionalBits**: The number of bits used to store the fractional part of coordinates in
   the fixed-point encoding.
6. **flags**: A bit field containing flags.
   - `0x1`: whether the splat was trained with [antialiasing](https://niujinshuchong.github.io/mip-splatting/).
7. **sh1Bits** (v3+): Number of quantization bits for SH degree 1 coefficients (1-8).
8. **shRestBits** (v3+): Number of quantization bits for SH degree 2+ coefficients (1-8).
9. **shMin** (v3+): Minimum SH coefficient value used for quantization scaling.
10. **shMax** (v3+): Maximum SH coefficient value used for quantization scaling.
11. **hasSafeOrbit** (v4+): Whether safe orbit camera data is present (0 or 1).
12. **safeOrbitElevationMin** (v4+): Minimum elevation for safe orbit in radians.
13. **safeOrbitElevationMax** (v4+): Maximum elevation for safe orbit in radians.
14. **safeOrbitRadiusMin** (v4+): Minimum radius for safe orbit.
15. **reserved**: Reserved for future use. Must be 0.

### Positions

Positions are represented as `(x, y, z)` coordinates, each as a 24-bit fixed point signed integer.
The number of fractional bits is determined by the `fractionalBits` field in the header.

### Scales

Scales are represented as `(x, y, z)` components, each represented as an 8-bit log-encoded integer.

### Rotation

Rotations are represented as the `(x, y, z)` components of the normalized rotation quaternion. The
`w` component can be derived from the others and is not stored. Each components is encoded as an
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

**Version 3 Quantization (current):**
- Uses adaptive min/max scaling based on the actual range of SH coefficients in the data
- Configurable quantization bits for degree 1 (default: 5 bits) and higher degrees (default: 4 bits)
- Stores quantization parameters in the file header for accurate reconstruction

**Legacy Quantization (versions 1-2):**
- Assumes SH coefficients are in the [-1, 1] range
- Fixed 5 bits of precision for degree 1 and 4 bits for degrees 2, 3, and 4
- Maintained for backward compatibility
