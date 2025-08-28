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

## API

```
bool saveSpz(
   const GaussianCloud &gaussians, PackOptions &options, std::vector<uint8_t> *output);
```

Converts a cloud of Gaussians in `.spz` format to a vector of bytes.

   - `gaussians`: The Gaussians to save
   - `options`: Flags that control the packing behavior.
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
   - `options`: Flags that control the packing behavior.
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

## File Format

The .spz format is a gzipped stream of data consisting of a 16-byte header followed by the
gaussian data. This data is organized by attribute in the following order: positions,
alphas, colors, scales, rotations, spherical harmonics.

### Header

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
2. **version**: Currently, the only valid versions are 2 and 3
3. **numPoints**: The number of gaussians
4. **shDegree**: The degree of spherical harmonics. This must be between 0 and 3 (inclusive).
5. **fractionalBits**: The number of bits used to store the fractional part of coordinates in
   the fixed-point encoding.
6. **flags**: A bit field containing flags.
   - `0x1`: whether the splat was trained with [antialiasing](https://niujinshuchong.github.io/mip-splatting/).
7. **reserved**: Reserved for future use. Must be 0.

### Positions

Positions are represented as `(x, y, z)` coordinates, each as a 24-bit fixed point signed integer.
The number of fractional bits is determined by the `fractionalBits` field in the header.

### Scales

Scales are represented as `(x, y, z)` components, each represented as an 8-bit log-encoded integer.

### Rotation

In version 3, rotations are represented as the smallest three components of the normalized rotation quaternion, for optimal rotation accuracy.
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
9 (for degree 1), 24 (for degree 2), or 45 (for degree 3) coefficients per gaussian.

The coefficients for a gaussian are organized such that the color channel is the inner (faster
varying) axis, and the coefficient is the outer (slower varying) axis, i.e. for degree 1,
the order of the 9 values is:

```
sh1n1_r, sh1n1_g, sh1n1_b, sh10_r, sh10_g, sh10_b, sh1p1_r, sh1p1_g, sh1p1_b
```

Each coefficient is represented as an 8-bit signed integer. Additional quantization can be performed
to attain a higher compression ratio. This library currently uses 5 bits of precision for degree 0
and 4 bits of precision for degrees 1 and 2, but this may be changed in the future without breaking
backwards compatibility.


## Python Bindings

The SPZ library provides Python bindings built with [nanobind](https://nanobind.readthedocs.io/) that offer a convenient interface for loading, manipulating, and saving 3D Gaussian splats from Python.

### Installation
```bash
git clone https://github.com/nianticlabs/spz.git
cd spz
pip install .
```

### API Overview

The Python API closely mirrors the C++ API with additional Python-friendly features:

- **Automatic type conversion**: Accepts numpy arrays with various dtypes (int32, float64, etc.) and converts them to float32
- **Memory safety**: All data is safely copied between Python and C++

### Basic Usage Examples

#### Loading and Saving SPZ Files

```python
import spz

# Load a .spz file
cloud = spz.load_spz("samples/hornedlizard.spz")
print(f"Loaded {cloud.num_points} gaussians with SH degree {cloud.sh_degree}")

# Save to a new file
options = spz.PackOptions()
options.from_coord = spz.CoordinateSystem.RUB
spz.save_spz(cloud, options, "output.spz")
```

#### Converting PLY to SPZ

```python
import spz

# Load a PLY file and convert to compressed SPZ format
# PLY files typically use RDF coordinates (right-handed, Y-down, Z-forward)
unpack_options = spz.UnpackOptions()
unpack_options.to_coord = spz.CoordinateSystem.RUB  # Convert to RUB for Three.js compatibility

# Load the PLY file
cloud = spz.load_splat_from_ply("input.ply", unpack_options)
print(f"Loaded {cloud.num_points} gaussians from PLY")

# Save as compressed SPZ format
pack_options = spz.PackOptions()
pack_options.from_coord = spz.CoordinateSystem.RUB  # Data is now in RUB coordinates
spz.save_spz(cloud, pack_options, "output.spz")

# Check compression ratio
import os
ply_size = os.path.getsize("input.ply")
spz_size = os.path.getsize("output.spz")
compression_ratio = ply_size / spz_size
print(f"Compression ratio: {compression_ratio:.1f}x smaller ({ply_size} → {spz_size} bytes)")
```

#### Creating and Manipulating Gaussian Clouds

```python
import spz
import numpy as np

# Create a new Gaussian cloud
cloud = spz.GaussianCloud()
cloud.sh_degree = 1
cloud.antialiased = True

# Choose point count for data you will assign
num_points = 100

# Set positions first (defines num_points)
positions = np.random.randn(num_points * 3).astype(np.float32)
cloud.positions = positions
assert cloud.num_points == num_points  # num_points is derived, read‑only

# Set scales (3 floats per point: log-scale factors)
scales = np.random.randn(num_points * 3).astype(np.float32)
cloud.scales = scales

# Set rotations (4 floats per point: quaternion x, y, z, w)
rotations = np.random.randn(num_points * 4).astype(np.float32)
cloud.rotations = rotations

# Set alphas (1 float per point: opacity before sigmoid)
alphas = np.random.randn(num_points).astype(np.float32)
cloud.alphas = alphas

# Set colors (3 floats per point: RGB)
colors = np.random.rand(num_points * 3).astype(np.float32)
cloud.colors = colors

# Set spherical harmonics (9 floats per point for degree 1)
sh_coeffs = np.random.randn(num_points * 9).astype(np.float32)
cloud.sh = sh_coeffs

# Calculate median volume
median_vol = cloud.median_volume()
print(f"Median volume: {median_vol}")

# Apply coordinate transformation
cloud.rotate_180_deg_about_x()  # Converts between RUB and RDF coordinates
```

#### Coordinate System Conversions

```python
import spz

# All available coordinate systems
print("Available coordinate systems:")
for coord_sys in [spz.CoordinateSystem.UNSPECIFIED, spz.CoordinateSystem.LDB, 
                  spz.CoordinateSystem.RDB, spz.CoordinateSystem.LUB, 
                  spz.CoordinateSystem.RUB, spz.CoordinateSystem.LDF,
                  spz.CoordinateSystem.RDF, spz.CoordinateSystem.LUF, 
                  spz.CoordinateSystem.RUF]:
    print(f"  {coord_sys}")

# Load PLY (typically RDF) and convert to Unity coordinates (RUF)
unpack_options = spz.UnpackOptions()
unpack_options.to_coord = spz.CoordinateSystem.RUF
cloud = spz.load_splat_from_ply("ply_file.ply", unpack_options)

# Save for Three.js (RUB coordinates)
pack_options = spz.PackOptions()
pack_options.from_coord = spz.CoordinateSystem.RUF  # Current data is in RUF
spz.save_spz(cloud, pack_options, "threejs_output.spz")  # Will be converted to RUB internally
```

#### In-place conversions on an existing GaussianCloud

```python
import spz
import numpy as np

cloud = spz.GaussianCloud()
cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
cloud.rotations = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)

# Convert from RUB to RDF (flips Y and Z axes)
cloud.convert_coordinates(spz.CoordinateSystem.RUB, spz.CoordinateSystem.RDF)

# Convert back to RUB
cloud.convert_coordinates(spz.CoordinateSystem.RDF, spz.CoordinateSystem.RUB)
```

#### Advanced Usage with NumPy Integration

```python
import spz
import numpy as np

# Load existing cloud
cloud = spz.load_spz("samples/hornedlizard.spz")

# Access data as NumPy arrays (always returns float32)
positions = cloud.positions  # Shape: (num_points * 3,)
scales = cloud.scales       # Shape: (num_points * 3,)
rotations = cloud.rotations # Shape: (num_points * 4,)
alphas = cloud.alphas       # Shape: (num_points,)
colors = cloud.colors       # Shape: (num_points * 3,)
sh_coeffs = cloud.sh        # Shape: (num_points * sh_coeffs_per_point,)

# Reshape for easier manipulation
positions_3d = positions.reshape(-1, 3)  # Shape: (num_points, 3)
colors_rgb = colors.reshape(-1, 3)       # Shape: (num_points, 3)

# Modify data
positions_3d[:, 2] += 1.0  # Move all points up by 1 unit in Z
colors_rgb[:, 0] *= 0.5    # Reduce red channel by half

# Update the cloud (automatic type conversion)
cloud.positions = positions_3d.flatten()
cloud.colors = colors_rgb.flatten()

# Save modified cloud
spz.save_spz(cloud, spz.PackOptions(), "modified.spz")
```

### Python API invariants and validations

The Python bindings enforce consistency across fields and provide clear errors:

- num_points is read‑only and derived from positions.size() / 3.
  - Set positions first to establish the point count.
- positions: length must be a multiple of 3. Setting positions updates num_points.
- scales: length must be a multiple of 3; if num_points > 0, length must equal num_points * 3.
- rotations: length must be a multiple of 4; if num_points > 0, length must equal num_points * 4.
- alphas: if num_points > 0, length must equal num_points.
- colors: length must be a multiple of 3; if num_points > 0, length must equal num_points * 3.
- sh_degree: must be in [0, 3]. Set sh_degree before assigning sh.
- sh:
  - If sh_degree == 0, sh must be empty.
  - Otherwise, length must be a multiple of (((sh_degree + 1)^2 − 1) * 3).
  - If num_points > 0, length must equal num_points * (((sh_degree + 1)^2 − 1) * 3).
- Dtypes: numeric arrays are accepted and converted to float32; non‑numeric arrays raise TypeError.
- All arrays must be C‑contiguous; non‑contiguous inputs will be copied by nanobind.

### Data Layout

The Python bindings maintain the same data layout as the C++ library:

- **Positions**: `[x1, y1, z1, x2, y2, z2, ...]`
- **Scales**: `[sx1, sy1, sz1, sx2, sy2, sz2, ...]` (log-scale)
- **Rotations**: `[x1, y1, z1, w1, x2, y2, z2, w2, ...]` (quaternions)
- **Alphas**: `[a1, a2, a3, ...]` (before sigmoid activation)
- **Colors**: `[r1, g1, b1, r2, g2, b2, ...]` (base RGB)
- **Spherical Harmonics**: Coefficient-major order, e.g., for degree 1:
  `[sh1n1_r, sh1n1_g, sh1n1_b, sh10_r, sh10_g, sh10_b, sh1p1_r, sh1p1_g, sh1p1_b, ...]`

### Type Safety

The Python bindings provide automatic type conversion while maintaining safety:

- ✅ **Accepts**: `int32`, `float64`, `uint8`, etc. → automatically converts to `float32`
- ❌ **Rejects**: `string`, `complex`, `object` arrays → raises `TypeError`
- 🔄 **Preserves**: `float32` arrays → no conversion needed

### Testing and Development

The Python bindings include a comprehensive test suite that covers all API functionality:

```bash

# Install pytest and scipy
pip install pytest scipy

# Run the test suite
python -m pytest tests/python/

### Requirements

- Python 3.8+
- NumPy (automatically installed)
- For development: pytest, scipy (for testing)
