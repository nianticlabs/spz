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
- `version`: Version of the packed format
- `sh1Bits`: Number of quantization bits for SH degree 1 coefficients (default: 5, range: 1-8)
- `shRestBits`: Number of quantization bits for SH degree 2+ coefficients (default: 4, range: 1-8)

## File Format

The .spz format is a gzipped stream of data consisting of a variable-size header followed by the
gaussian data. This data is organized by attribute in the following order: positions,
alphas, colors, scales, rotations, spherical harmonics.

### Header

**Version 2 (current):**
```c
struct PackedGaussiansHeader {
  uint32_t magic;
  uint32_t version;
  uint32_t numPoints;
  uint8_t shDegree;
  uint8_t fractionalBits;
  uint8_t flags;
  uint8_t v2Padding;
};
```

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

**Quantization:**

By default, a fixed-precision quantization is adopted, where we
- Assumes SH coefficients are in the [-1, 1] range
- Fixed 5 bits of precision for degree 1 and 4 bits for degrees 2, 3, and 4
- Maintained for backward compatibility

### Extensions

### Spherical Harmonics Quantization

With extension `SPZ_ADOBE_sh_quantization`, SPZ supports configurable spherical harmonics quantization with the following improvements:

1. **Adaptive Range Scaling**: Instead of assuming SH coefficients are in the [-1, 1] range, the format now computes the actual min/max values from the data and uses the full quantization range for optimal precision.

2. **Configurable Bit Precision**: The number of quantization bits for SH degree 1 and higher degree coefficients can be configured (1-8 bits), allowing users to trade off between file size and quality. Default bits are 5 for SH degree 1 and 4 for higher ones.

3. **Backward Compatibility**: The library maintains full backward compatibility when this extension is unused.

**Attributes**

This extension has the following attributes and default values:

```
  uint8_t sh1Bits = 5;     // Bits for SH degree 1 coefficients
  uint8_t shRestBits = 4;  // Bits for SH degree 2+ coefficients
  float shMin = -1.0f;     // Minimum SH coefficient value used for quantization
  float shMax = 1.0f;      // Maximum SH coefficient value used for quantization
```

### Camera Orbit Limitation

With extension `SPZ_ADOBE_safe_orbit_camera`, SPZ supports storing camera limits which can be used to restrict the view in a render. This extension includes:

**Attributes**

This extension has the following attributes and default values:

```
  float safeOrbitElevationMin = 0.0f;  // Minimum elevation for safe orbit (radians)
  float safeOrbitElevationMax = 0.0f;  // Maximum elevation for safe orbit (radians)
  float safeOrbitRadiusMin = 0.0f;     // Minimum radius for safe orbit
```
