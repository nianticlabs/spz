## Python Bindings

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
