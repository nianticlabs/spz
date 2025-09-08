import os
import tempfile

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

import spz

# -----------------------------------------------------------------------------
# Helper functions using SciPy for quaternion math.
# -----------------------------------------------------------------------------

def normalized(v):
    """Return the normalized version of v (as a numpy array)."""
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-8:
        return v
    return v / n

def axis_angle_quat(angle_axis):
    """
    Convert an axis–angle vector (where the angle is the norm) into a quaternion.
    SciPy’s Rotation.from_rotvec returns a quaternion in [x, y, z, w] order.
    We convert that to [w, x, y, z] order.
    """
    r = Rotation.from_rotvec(angle_axis)
    q = r.as_quat()  # [x, y, z, w]
    return np.concatenate(([q[3]], q[:3]))

def times(a, b):
    """
    Overloaded multiplication:
      - If b is a scalar, return a * b.
      - If a and b are both 4-element arrays, treat them as quaternions (in [w,x,y,z] order)
        and return their product (also in [w,x,y,z] order).
      - If a is a quaternion (4-element array in [w,x,y,z] order) and b is a 3-vector,
        rotate b by the quaternion.
      - Otherwise, perform element–wise multiplication.
    """
    a = np.array(a, dtype=float)
    # Scalar multiplication:
    if isinstance(b, (int, float)):
        return a * b

    b_arr = np.array(b, dtype=float)
    # Quaternion multiplication: both are 4-element arrays.
    if a.shape == (4,) and b_arr.shape == (4,):
        r1 = Rotation.from_quat([a[1], a[2], a[3], a[0]])
        r2 = Rotation.from_quat([b_arr[1], b_arr[2], b_arr[3], b_arr[0]])
        r3 = r1 * r2
        q = r3.as_quat()  # SciPy returns [x, y, z, w]
        return np.array([q[3], q[0], q[1], q[2]])
    # If a is a quaternion and b is a 3-vector, rotate b.
    if a.shape == (4,) and b_arr.shape == (3,):
        r = Rotation.from_quat([a[1], a[2], a[3], a[0]])
        return r.apply(b_arr)
    # Otherwise, element–wise multiplication.
    return a * b_arr

# -----------------------------------------------------------------------------
# Other helper functions for the tests.
# -----------------------------------------------------------------------------


def read_file(path):
    """Read the entire file as a string (binary read, then decode)."""
    with open(path, "rb") as f:
        return f.read().decode("utf-8", errors="replace")

def make_test_gaussian_cloud(include_sh):
    """
    Create a GaussianCloud with two splats for testing.
    If include_sh is True then spherical harmonics (SH) are added and sh_degree is set to 3.
    """
    # cloud = spz.GaussianCloud(
    #     num_points=2,
    #     antialiased=True,
    #     positions=[0, 0.1, -0.2, 0.3, 0.4, 0.5],
    #     scales=[-3, -2, -1.5, -1, 0, 0.1],
    #     rotations=[-0.5, 0.2, 1, -0.2, 0.1, -0.4, -0.3, 0.5],
    #     alphas=[-1.0, 1.0],
    #     colors=[-1, 0, 1, -0.5, 0.5, 0.1],
    # )
    cloud = spz.GaussianCloud()
    cloud.antialiased = True
    cloud.positions = np.array([0, 0.1, -0.2, 0.3, 0.4, 0.5], dtype=float)
    cloud.scales = np.array([-3, -2, -1.5, -1, 0, 0.1], dtype=float)
    cloud.rotations = np.array([-0.5, 0.2, 1, -0.2, 0.1, -0.4, -0.3, 0.5], dtype=float)
    cloud.alphas = np.array([-1.0, 1.0], dtype=float)
    cloud.colors = np.array([-1, 0, 1, -0.5, 0.5, 0.1], dtype=float)
    if include_sh:
        # Degree 3 -> 45 coeffs per point × 2 points = 90
        cloud.sh_degree = 3
        cloud.sh = np.array([i / 45.0 - 1.0 for i in range(90)], dtype=float)
    else:
        cloud.sh_degree = 0
        cloud.sh = np.array([], dtype=float)
    return cloud

# -----------------------------------------------------------------------------
# Epsilon constants (from the C++ tests)
# -----------------------------------------------------------------------------

SH_4BIT_EPSILON = 2.0 / 32.0 + 0.5 / 255.0
SH_5BIT_EPSILON = 2.0 / 64.0 + 0.5 / 255.0

# -----------------------------------------------------------------------------
# Test functions
# -----------------------------------------------------------------------------

def test_save_load_packed_format():
    """Test saving and loading SPZ format with compression and precision checks."""
    src = make_test_gaussian_cloud(include_sh=True)
    filename = os.path.join(tempfile.gettempdir(), "SplatIOTest_SaveLoad.spz")
    assert spz.save_spz(src, spz.PackOptions(), filename) is True

    dst = spz.load_spz(filename, spz.UnpackOptions())
    assert dst.num_points == 2
    assert dst.sh_degree == 3
    assert dst.antialiased is True

    # Compare positions and scales.
    np.testing.assert_allclose(dst.positions, src.positions, atol=1 / 2048.0)
    np.testing.assert_allclose(dst.scales, src.scales, atol=1 / 32.0)

    # Check rotations: extract the first two quaternions (each 4 numbers) and normalize.
    q0 = np.array(dst.rotations[0:4], dtype=float)
    q1 = np.array(dst.rotations[4:8], dtype=float)
    orig_q0 = normalized(np.array(src.rotations[0:4], dtype=float))
    orig_q1 = normalized(np.array(src.rotations[4:8], dtype=float))

    assert np.isclose(np.linalg.norm(q0), 1.0, atol=1e-6)
    assert np.isclose(np.linalg.norm(q1), 1.0, atol=1e-6)

    v1 = np.array([3.0, -2.0, 0.2])
    v2 = np.array([-1.0, 0.5, -3.0])
    for q, orig_q in [(q0, orig_q0), (q1, orig_q1)]:
        for v in [v1, v2]:
            a = times(q, v)
            b = times(orig_q, v)
            cosine = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
            assert np.isclose(cosine, 1.0, atol=1e-4)

    np.testing.assert_allclose(dst.alphas, src.alphas, atol=0.01)
    np.testing.assert_allclose(dst.sh, src.sh, atol=SH_4BIT_EPSILON)
    # Check degree‑1 SH (first 9 coefficients) with extra precision.
    np.testing.assert_allclose(dst.sh[0:9], src.sh[0:9], atol=SH_5BIT_EPSILON)
    np.testing.assert_allclose(dst.sh[45:45 + 9], src.sh[45:45 + 9], atol=SH_5BIT_EPSILON)

def test_save_load_packed_format_large_splat():
    """Test saving and loading large SPZ files with many points."""
    num_points = 50000
    # src = spz.GaussianCloud(numPoints=num_points, shDegree=3)
    src = spz.GaussianCloud()
    src.sh_degree = 3

    rng = np.random.default_rng(1)
    src.positions = (rng.uniform(0.0, 1.0, size=(num_points, 3)) * 2.0 - 1.0).flatten()
    src.scales = (rng.uniform(0.0, 1.0, size=(num_points, 3)) - 1.0).flatten()
    src.rotations = (rng.uniform(0.0, 1.0, size=(num_points, 4)) * 2.0 - 1.0).flatten()
    src.colors = rng.uniform(0.0, 1.0, size=(num_points, 3)).flatten()
    src.alphas = rng.uniform(0.0, 1.0, size=num_points)
    src.sh = (rng.uniform(0.0, 1.0, size=(num_points, 45)) - 0.5).flatten()

    filename = os.path.join(tempfile.gettempdir(), "large_splat.spz")
    assert spz.save_spz(src, spz.PackOptions(), filename) is True

    dst = spz.load_spz(filename, spz.UnpackOptions())
    assert dst.num_points == src.num_points
    assert dst.sh_degree == src.sh_degree
    np.testing.assert_allclose(dst.positions, src.positions, atol=1 / 2048.0)
    np.testing.assert_allclose(dst.scales, src.scales, atol=1 / 16.0)
    assert len(dst.rotations) == len(src.rotations)
    np.testing.assert_allclose(dst.alphas, src.alphas, atol=0.01)
    sh_epsilon = 2.0 / 32.0 + 1.0 / 255.0
    np.testing.assert_allclose(dst.sh, src.sh, atol=sh_epsilon)

def test_sh_encoding_for_zeros_and_edges():
    """Test spherical harmonics encoding for edge values and zeros."""
    # src = spz.GaussianCloud(
    #     numPoints=1,
    #     shDegree=1,
    #     positions=np.zeros(3),
    #     scales=np.zeros(3),
    #     rotations=np.array([0, 0, 0, 1]),
    #     alphas=np.array([0]),
    #     colors=np.zeros(3),
    #     sh=np.array([-0.01, 0.0, 0.01, -1.0, -0.99, -0.95, 0.95, 0.99, 1.0]),
    # )
    src = spz.GaussianCloud()
    src.sh_degree = 1
    src.positions = np.zeros(3, dtype=float)
    src.scales = np.zeros(3, dtype=float)
    src.rotations = np.array([0, 0, 0, 1], dtype=float)
    src.alphas = np.array([0.0], dtype=float)
    src.colors = np.zeros(3, dtype=float)
    src.sh = np.array([-0.01, 0.0, 0.01, -1.0, -0.99, -0.95, 0.95, 0.99, 1.0], dtype=float)

    filename = os.path.join(tempfile.gettempdir(), "test_sh_encoding_for_zeros_and_edges.spz")
    assert spz.save_spz(src, spz.PackOptions(), filename) is True
    dst = spz.load_spz(filename, spz.UnpackOptions())
    assert dst.num_points == 1
    assert dst.sh_degree == 1
    expected_sh = np.array([0.0, 0.0, 0.0, -1.0, -1.0, -0.9375, 0.9375, 0.9922, 0.9922])
    np.testing.assert_allclose(dst.sh, expected_sh, atol=2e-5)

@pytest.mark.parametrize("include_sh", [False, True])
def test_save_load_ply(include_sh):
    """Test saving and loading PLY files with and without spherical harmonics."""
    src = make_test_gaussian_cloud(include_sh=include_sh)
    filename = os.path.join(tempfile.gettempdir(), "SplatIOTest_SaveLoad.ply")
    assert spz.save_splat_to_ply(src, spz.PackOptions(), filename) is True

    ply = read_file(filename)
    expected_header = "ply\nformat binary_little_endian 1.0\nelement vertex 2\n"
    assert ply.startswith(expected_header)

    dst = spz.load_splat_from_ply(filename, spz.UnpackOptions())
    assert dst.num_points == 2
    assert dst.sh_degree == (3 if include_sh else 0)
    assert np.array_equal(dst.positions, src.positions)
    assert np.array_equal(dst.scales, src.scales)
    assert np.array_equal(dst.rotations, src.rotations)
    assert np.array_equal(dst.alphas, src.alphas)
    assert np.array_equal(dst.colors, src.colors)
    if include_sh:
        assert np.array_equal(dst.sh, src.sh)
    else:
        assert len(dst.sh) == 0

def test_coordinate_system_enum():
    """Test that all coordinate system enum values are available and unique."""
    # Test all enum values exist
    assert hasattr(spz, 'CoordinateSystem')
    assert hasattr(spz, 'UNSPECIFIED')
    assert hasattr(spz, 'LDB')
    assert hasattr(spz, 'RDB')
    assert hasattr(spz, 'LUB')
    assert hasattr(spz, 'RUB')
    assert hasattr(spz, 'LDF')
    assert hasattr(spz, 'RDF')
    assert hasattr(spz, 'LUF')
    assert hasattr(spz, 'RUF')
    
    # Test that all enum values are unique
    enum_values = [
        spz.UNSPECIFIED, spz.LDB, spz.RDB, spz.LUB, spz.RUB,
        spz.LDF, spz.RDF, spz.LUF, spz.RUF
    ]
    assert len(enum_values) == len(set(enum_values))
    
    # Test that coordinate systems can be used in options
    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.LDB
    assert pack_opts.from_coord == spz.LDB
    
    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RUF
    assert unpack_opts.to_coord == spz.RUF


def test_pack_options_mutability():
    """Test that PackOptions can be created and modified."""
    opts = spz.PackOptions()
    
    # Test default initialization
    assert opts.from_coord == spz.UNSPECIFIED
    
    # Test setting different coordinate systems
    opts.from_coord = spz.LDB
    assert opts.from_coord == spz.LDB
    
    opts.from_coord = spz.RUF
    assert opts.from_coord == spz.RUF


def test_unpack_options_mutability():
    """Test that UnpackOptions can be created and modified."""
    opts = spz.UnpackOptions()
    
    # Test default initialization
    assert opts.to_coord == spz.UNSPECIFIED
    
    # Test setting different coordinate systems
    opts.to_coord = spz.RDB
    assert opts.to_coord == spz.RDB
    
    opts.to_coord = spz.LUF
    assert opts.to_coord == spz.LUF


def test_gaussian_cloud_initialization():
    """Test GaussianCloud initialization and basic properties."""
    cloud = spz.GaussianCloud()
    
    # Test default values
    assert cloud.num_points == 0
    assert cloud.sh_degree == 0
    assert cloud.antialiased == False
    
    # Test that arrays are initially empty
    assert len(cloud.positions) == 0
    assert len(cloud.scales) == 0
    assert len(cloud.rotations) == 0
    assert len(cloud.alphas) == 0
    assert len(cloud.colors) == 0
    assert len(cloud.sh) == 0


def test_gaussian_cloud_property_setting():
    """Test setting properties on GaussianCloud."""
    cloud = spz.GaussianCloud()
    
    # num_points is read-only; attempting to set should fail
    with pytest.raises(AttributeError):
        cloud.num_points = 5
    # Test setting scalar properties
    cloud.sh_degree = 2
    cloud.antialiased = True
    
    assert cloud.sh_degree == 2
    assert cloud.antialiased == True
    
    # Test setting array properties
    positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.positions = positions
    np.testing.assert_array_equal(cloud.positions, positions)
    assert cloud.num_points == 1
    
    scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.scales = scales
    np.testing.assert_array_equal(cloud.scales, scales)
    
    rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.rotations = rotations
    np.testing.assert_array_equal(cloud.rotations, rotations)
    
    alphas = np.array([0.5], dtype=np.float32)
    cloud.alphas = alphas
    np.testing.assert_array_equal(cloud.alphas, alphas)
    
    colors = np.array([1.0, 0.0, 0.0], dtype=np.float32)
    cloud.colors = colors
    np.testing.assert_array_equal(cloud.colors, colors)
    
    # For sh_degree=2, per-point SH length must be 24
    sh = np.zeros(24, dtype=np.float32)
    cloud.sh = sh
    np.testing.assert_array_equal(cloud.sh, sh)


def test_gaussian_cloud_array_dtype_handling():
    """Test that GaussianCloud properly handles different array dtypes."""
    cloud = spz.GaussianCloud()
    
    # Test with different numpy dtypes - should convert to float32
    for dtype in [np.float64, np.int32, np.float32]:
        data = np.array([1.0, 2.0, 3.0], dtype=dtype)
        cloud.positions = data
        # The returned array should always be float32
        assert cloud.positions.dtype == np.float32
        np.testing.assert_array_equal(cloud.positions, [1.0, 2.0, 3.0])
    
    # Test that non-numeric types are rejected
    with pytest.raises(TypeError, match="incompatible function arguments"):
        cloud.positions = np.array(['a', 'b', 'c'], dtype=np.str_)
    
    # Test that complex types are rejected
    with pytest.raises(TypeError, match="incompatible function arguments"):
        cloud.positions = np.array([1+2j, 3+4j, 5+6j], dtype=np.complex64)


def test_gaussian_cloud_rotate_180_deg_about_x():
    """Test the rotate_180_deg_about_x method."""
    cloud = spz.GaussianCloud()
    
    # Set up a test point with known coordinates
    original_pos = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    original_rot = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)  # Non-identity quaternion
    
    cloud.positions = original_pos.copy()
    cloud.rotations = original_rot.copy()
    
    # Apply rotation - this converts from RUB to RDF coordinates
    cloud.rotate_180_deg_about_x()
    
    # According to the C++ implementation, this converts between RUB and RDF
    # RUB: Right Up Back, RDF: Right Down Front
    # This should flip the Y and Z coordinates and corresponding quaternion components
    
    # Check that positions have been transformed (Y and Z should be flipped)
    expected_pos = np.array([1.0, -2.0, -3.0], dtype=np.float32)
    np.testing.assert_array_almost_equal(cloud.positions, expected_pos, decimal=5)
    
    # Check that rotations have been transformed
    # The quaternion should have Y and Z components flipped according to coordinate conversion
    expected_rot = np.array([0.1, -0.2, -0.3, 0.9], dtype=np.float32)
    np.testing.assert_array_almost_equal(cloud.rotations, expected_rot, decimal=5)


def test_gaussian_cloud_median_volume():
    """Test the median_volume method calculates volume correctly."""
    cloud = spz.GaussianCloud()
    
    # Test empty cloud - should return 0.01 according to C++ implementation
    np.testing.assert_almost_equal(cloud.median_volume(), 0.01, decimal=5)
    
    # Test with actual points (3 points)
    cloud.positions = np.zeros(3 * 3, dtype=np.float32)
    
    # Set up scales for volume calculation (log scale)
    # Volume = 4/3 * pi * exp(scale_sum)
    # Use known values: scales [-1,-1,-1], [0,0,0], [1,1,1]
    # Scale sums: -3, 0, 3
    # Median scale sum: 0
    # Expected median volume: 4/3 * pi * exp(0) = 4/3 * pi
    cloud.scales = np.array([
        -1.0, -1.0, -1.0,  # First gaussian: scale sum = -3
        0.0, 0.0, 0.0,     # Second gaussian: scale sum = 0
        1.0, 1.0, 1.0      # Third gaussian: scale sum = 3
    ], dtype=np.float32)
    
    median_vol = cloud.median_volume()
    expected_vol = (4.0 / 3.0) * np.pi * np.exp(0.0)  # exp(0) = 1
    np.testing.assert_almost_equal(median_vol, expected_vol, decimal=5)
    
    # Test with 5 points to verify median calculation
    cloud.positions = np.zeros(5 * 3, dtype=np.float32)
    cloud.scales = np.array([
        -2.0, -2.0, -2.0,  # scale sum = -6
        -1.0, -1.0, -1.0,  # scale sum = -3
        0.0, 0.0, 0.0,     # scale sum = 0 (median)
        1.0, 1.0, 1.0,     # scale sum = 3
        2.0, 2.0, 2.0      # scale sum = 6
    ], dtype=np.float32)
    
    median_vol = cloud.median_volume()
    expected_vol = (4.0 / 3.0) * np.pi * np.exp(0.0)  # median scale sum is 0
    np.testing.assert_almost_equal(median_vol, expected_vol, decimal=5)


def test_coordinate_system_conversion():
    """Test coordinate system conversion during file I/O with actual verification."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1
    cloud.antialiased = False
    
    # Create test data in RUB coordinates (Right Up Back)
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.rotations = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.sh = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], dtype=np.float32)
    
    filename = os.path.join(tempfile.gettempdir(), "coord_conversion_test.spz")
    
    # Save as RUB (no conversion from RUB to RUB)
    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RUB
    assert spz.save_spz(cloud, pack_opts, filename) is True
    
    # Load with conversion to RDF (Right Down Front)
    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RDF
    loaded_cloud = spz.load_spz(filename, unpack_opts)
    
    # Verify the cloud was loaded successfully
    assert loaded_cloud.num_points == cloud.num_points
    assert loaded_cloud.sh_degree == cloud.sh_degree
    assert loaded_cloud.antialiased == cloud.antialiased
    
    # RUB to RDF conversion should flip Y and Z coordinates
    # According to the C++ implementation, this affects positions and rotations
    expected_pos = np.array([1.0, -2.0, -3.0], dtype=np.float32)
    np.testing.assert_allclose(loaded_cloud.positions, expected_pos, atol=1/2048.0)
    
    # Quaternion Y and Z components should be flipped (but not W)
    expected_rot = np.array([0.1, -0.2, -0.3, 0.9], dtype=np.float32)
    # Normalize both quaternions for comparison since they get normalized during packing
    loaded_rot_norm = loaded_cloud.rotations / np.linalg.norm(loaded_cloud.rotations)
    expected_rot_norm = expected_rot / np.linalg.norm(expected_rot)
    np.testing.assert_allclose(loaded_rot_norm, expected_rot_norm, atol=1e-3)
    
    # Test conversion from RDF to LUF (more complex transformation)
    pack_opts.from_coord = spz.RDF
    unpack_opts.to_coord = spz.LUF
    
    # Create new test cloud in RDF coordinates
    cloud2 = spz.GaussianCloud()
    cloud2.sh_degree = 0
    cloud2.antialiased = False
    cloud2.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud2.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud2.rotations = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)
    cloud2.alphas = np.array([0.5], dtype=np.float32)
    cloud2.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud2.sh = np.array([], dtype=np.float32)
    
    filename2 = os.path.join(tempfile.gettempdir(), "coord_conversion_test2.spz")
    assert spz.save_spz(cloud2, pack_opts, filename2) is True
    
    loaded_cloud2 = spz.load_spz(filename2, unpack_opts)
    assert loaded_cloud2.num_points == cloud2.num_points
    assert loaded_cloud2.sh_degree == cloud2.sh_degree
    
    # The transformation should affect the coordinates
    # RDF (Right Down Front) to LUF (Left Up Front) should flip X and Y
    expected_pos2 = np.array([-1.0, -2.0, 3.0], dtype=np.float32)
    np.testing.assert_allclose(loaded_cloud2.positions, expected_pos2, atol=1/2048.0)


def test_spherical_harmonics_degree_coefficients():
    """Test that spherical harmonics coefficients match expected counts for different degrees."""
    cloud = spz.GaussianCloud()
    
    # Test degree 0 (no SH coefficients)
    cloud.sh_degree = 0
    cloud.sh = np.array([], dtype=np.float32)
    expected_coeffs = 0
    assert len(cloud.sh) == expected_coeffs
    
    # Test degree 1 (9 coefficients per point: 3 coeffs × 3 channels)
    cloud.sh_degree = 1
    cloud.sh = np.array([0.0] * 9, dtype=np.float32)
    expected_coeffs = 9
    assert len(cloud.sh) == expected_coeffs
    
    # Test degree 2 (24 coefficients per point: 8 coeffs × 3 channels)
    cloud.sh_degree = 2
    cloud.sh = np.array([0.0] * 24, dtype=np.float32)
    expected_coeffs = 24
    assert len(cloud.sh) == expected_coeffs
    
    # Test degree 3 (45 coefficients per point: 15 coeffs × 3 channels)
    cloud.sh_degree = 3
    cloud.sh = np.array([0.0] * 45, dtype=np.float32)
    expected_coeffs = 45
    assert len(cloud.sh) == expected_coeffs
    
    # Test multiple points
    cloud.positions = np.zeros(2 * 3, dtype=np.float32)
    cloud.sh_degree = 1
    cloud.sh = np.array([0.0] * 18, dtype=np.float32)  # 2 points × 9 coeffs
    assert len(cloud.sh) == 18


def test_spherical_harmonics_coordinate_transformation():
    """Test that spherical harmonics are properly transformed during coordinate conversion."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1  # 9 coefficients per point
    cloud.antialiased = False
    
    # Set up test data with non-zero SH coefficients
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    
    # Set up SH coefficients for degree 1 (3 coefficients × 3 channels)
    # The ordering is: sh1n1_r, sh1n1_g, sh1n1_b, sh10_r, sh10_g, sh10_b, sh1p1_r, sh1p1_g, sh1p1_b
    original_sh = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], dtype=np.float32)
    cloud.sh = original_sh.copy()
    
    filename = os.path.join(tempfile.gettempdir(), "sh_coord_test.spz")
    
    # Save as RUB and load as RDF (180 degree rotation about X)
    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RUB
    assert spz.save_spz(cloud, pack_opts, filename) is True
    
    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RDF
    loaded_cloud = spz.load_spz(filename, unpack_opts)
    
    # Verify basic properties
    assert loaded_cloud.num_points == 1
    assert loaded_cloud.sh_degree == 1
    assert len(loaded_cloud.sh) == 9
    
    # According to C++ implementation, coordinate conversion affects SH coefficients
    # The conversion flips certain coefficients based on the coordinate system transformation
    # This is a complex transformation, but we can verify that it's been applied
    
    # The SH coefficients should be different from the original (transformed)
    # Some coefficients should be flipped according to the coordinate system change
    assert not np.array_equal(loaded_cloud.sh, original_sh)
    
    # Verify the transformation is consistent - save and load again should give same result
    filename2 = os.path.join(tempfile.gettempdir(), "sh_coord_test2.spz")
    pack_opts2 = spz.PackOptions()
    pack_opts2.from_coord = spz.RDF
    assert spz.save_spz(loaded_cloud, pack_opts2, filename2) is True
    
    unpack_opts2 = spz.UnpackOptions()
    unpack_opts2.to_coord = spz.RDF
    loaded_cloud2 = spz.load_spz(filename2, unpack_opts2)
    
    np.testing.assert_array_almost_equal(loaded_cloud.sh, loaded_cloud2.sh, decimal=4)


def test_quaternion_normalization_during_packing():
    """Test that quaternions are normalized during SPZ packing/unpacking."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.antialiased = False
    
    # Set up test data with non-normalized quaternions
    cloud.positions = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], dtype=np.float32)
    cloud.alphas = np.array([0.5, 0.7], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], dtype=np.float32)
    cloud.sh = np.array([], dtype=np.float32)
    
    # Use non-normalized quaternions
    non_normalized_quats = np.array([
        2.0, 3.0, 4.0, 5.0,  # First quaternion (not normalized)
        1.0, 1.0, 1.0, 1.0   # Second quaternion (not normalized)
    ], dtype=np.float32)
    cloud.rotations = non_normalized_quats.copy()
    
    # Save and load the cloud
    filename = os.path.join(tempfile.gettempdir(), "quat_normalization_test.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), filename) is True
    
    loaded_cloud = spz.load_spz(filename, spz.UnpackOptions())
    
    # Verify that the quaternions are now normalized
    assert len(loaded_cloud.rotations) == 8  # 2 quaternions × 4 components
    
    # Check first quaternion is normalized
    q1 = loaded_cloud.rotations[0:4]
    q1_norm = np.linalg.norm(q1)
    assert abs(q1_norm - 1.0) < 1e-4
    
    # Check second quaternion is normalized
    q2 = loaded_cloud.rotations[4:8]
    q2_norm = np.linalg.norm(q2)
    assert abs(q2_norm - 1.0) < 1e-4
    
    # Verify the orientation is preserved (normalized quaternions should represent same rotation)
    # Calculate expected normalized quaternions
    expected_q1 = non_normalized_quats[0:4] / np.linalg.norm(non_normalized_quats[0:4])
    expected_q2 = non_normalized_quats[4:8] / np.linalg.norm(non_normalized_quats[4:8])
    
    # Due to compression, we need to allow for some tolerance
    # Also, quaternions can be negated and still represent the same rotation
    def quaternions_equivalent(q1, q2, tolerance=1e-2):
        return (np.allclose(q1, q2, atol=tolerance) or 
                np.allclose(q1, -q2, atol=tolerance))
    
    assert quaternions_equivalent(q1, expected_q1) or quaternions_equivalent(q1, -expected_q1)
    assert quaternions_equivalent(q2, expected_q2) or quaternions_equivalent(q2, -expected_q2)


def test_convert_coordinates_method():
    """Directly test in-place coordinate conversion binding."""
    cloud = spz.GaussianCloud()
    # Start with RUB data
    cloud.sh_degree = 0
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.rotations = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)  # x, y, z, w

    # Convert to RDF: should flip Y and Z for positions and quaternion components
    cloud.convert_coordinates(spz.RUB, spz.RDF)
    np.testing.assert_allclose(cloud.positions, np.array([1.0, -2.0, -3.0], dtype=np.float32), atol=1e-6)
    np.testing.assert_allclose(cloud.rotations, np.array([0.1, -0.2, -0.3, 0.9], dtype=np.float32), atol=1e-6)

    # Convert back to RUB: should restore original values
    cloud.convert_coordinates(spz.RDF, spz.RUB)
    np.testing.assert_allclose(cloud.positions, np.array([1.0, 2.0, 3.0], dtype=np.float32), atol=1e-6)
    np.testing.assert_allclose(cloud.rotations, np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32), atol=1e-6)

def test_edge_cases_empty_arrays():
    """Test handling of empty arrays and edge cases."""
    cloud = spz.GaussianCloud()
    
    # Test empty arrays
    empty_array = np.array([], dtype=np.float32)
    cloud.positions = empty_array
    cloud.scales = empty_array
    cloud.rotations = empty_array
    cloud.alphas = empty_array
    cloud.colors = empty_array
    cloud.sh = empty_array
    
    assert len(cloud.positions) == 0
    assert len(cloud.scales) == 0
    assert len(cloud.rotations) == 0
    assert len(cloud.alphas) == 0
    assert len(cloud.colors) == 0
    assert len(cloud.sh) == 0


def test_compression_precision_validation():
    """Test that compression maintains expected precision levels."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1
    cloud.antialiased = False
    
    # Test with specific values that test the compression boundaries
    # Positions: 12-bit fractional precision (1/4096 resolution)
    cloud.positions = np.array([1.0, -1.0, 0.5], dtype=np.float32)
    
    # Scales: 5-bit precision (1/32 resolution)
    cloud.scales = np.array([1.0, -1.0, 0.5], dtype=np.float32)
    
    # Rotations: quaternion with smallest-three compression
    cloud.rotations = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)
    
    # Alpha: 8-bit precision
    cloud.alphas = np.array([0.5], dtype=np.float32)
    
    # Colors: 8-bit precision
    cloud.colors = np.array([0.5, -0.5, 0.25], dtype=np.float32)
    
    # SH: 4-bit precision for most coefficients, 5-bit for degree-1
    cloud.sh = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], dtype=np.float32)
    
    filename = os.path.join(tempfile.gettempdir(), "compression_precision_test.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), filename) is True
    
    loaded_cloud = spz.load_spz(filename, spz.UnpackOptions())
    
    # Verify precision according to C++ implementation
    # Positions: 12-bit fractional (1/4096 ≈ 0.000244)
    np.testing.assert_allclose(loaded_cloud.positions, cloud.positions, atol=1/2048.0)
    
    # Scales: 5-bit precision (1/32 = 0.03125)
    np.testing.assert_allclose(loaded_cloud.scales, cloud.scales, atol=1/32.0)
    
    # Quaternions should be normalized and preserve orientation
    loaded_quat = loaded_cloud.rotations
    assert abs(np.linalg.norm(loaded_quat) - 1.0) < 1e-4
    
    # Alpha: 8-bit precision
    np.testing.assert_allclose(loaded_cloud.alphas, cloud.alphas, atol=0.01)
    
    # Colors: 8-bit precision
    np.testing.assert_allclose(loaded_cloud.colors, cloud.colors, atol=0.01)
    
    # SH: degree-1 coefficients (first 9) have 5-bit precision
    np.testing.assert_allclose(loaded_cloud.sh[0:9], cloud.sh[0:9], atol=SH_5BIT_EPSILON)
    
    # Verify the constants match the C++ implementation
    assert SH_4BIT_EPSILON == 2.0 / 32.0 + 0.5 / 255.0
    assert SH_5BIT_EPSILON == 2.0 / 64.0 + 0.5 / 255.0


def test_edge_cases_empty_cloud():
    """Test operations on empty GaussianCloud."""
    cloud = spz.GaussianCloud()
    
    # Test median_volume on empty cloud - should return 0.01 according to C++ implementation
    np.testing.assert_almost_equal(cloud.median_volume(), 0.01, decimal=5)
    
    # Test rotate_180_deg_about_x on empty cloud
    cloud.rotate_180_deg_about_x()  # Should not crash
    
    # Test saving empty cloud
    filename = os.path.join(tempfile.gettempdir(), "empty_cloud.spz")
    result = spz.save_spz(cloud, spz.PackOptions(), filename)
    assert result is True
    
    # Test loading empty cloud
    loaded_cloud = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded_cloud.num_points == 0
    assert loaded_cloud.sh_degree == 0
    assert len(loaded_cloud.positions) == 0


def test_performance_large_cloud():
    """Test performance with a large cloud and verify timing."""
    import time
    
    num_points = 10000
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 2
    
    # Create large arrays
    rng = np.random.default_rng(42)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 24).astype(np.float32)
    
    # Time the save operation
    filename = os.path.join(tempfile.gettempdir(), "large_cloud_performance.spz")
    start_time = time.time()
    result = spz.save_spz(cloud, spz.PackOptions(), filename)
    save_time = time.time() - start_time
    
    assert result is True
    assert save_time < 5.0  # Should complete within 5 seconds
    
    # Time the load operation
    start_time = time.time()
    loaded_cloud = spz.load_spz(filename, spz.UnpackOptions())
    load_time = time.time() - start_time
    
    assert loaded_cloud.num_points == num_points
    assert load_time < 5.0  # Should complete within 5 seconds


def test_error_handling_invalid_arrays():
    """Test error handling for invalid array inputs."""
    cloud = spz.GaussianCloud()
    
    # Test with non-float32 arrays that can't be converted
    with pytest.raises(TypeError):
        # Create a complex array that can't be converted to float32
        complex_array = np.array([1+2j, 3+4j], dtype=np.complex64)
        cloud.positions = complex_array

    # Shape validations for properties
    cloud.sh_degree = 2
    cloud.positions = np.zeros(3, dtype=np.float32)  # one point
    with pytest.raises(ValueError, match="positions length must be a multiple of 3"):
        cloud.positions = np.array([1, 2], dtype=np.float32)
    with pytest.raises(ValueError, match=r"scales length must equal num_points \* 3"):
        cloud.scales = np.zeros(6, dtype=np.float32)  # expects 3
    with pytest.raises(ValueError, match=r"rotations length must equal num_points \* 4"):
        cloud.rotations = np.zeros(8, dtype=np.float32)  # expects 4
    with pytest.raises(ValueError, match=r"colors length must equal num_points \* 3"):
        cloud.colors = np.zeros(6, dtype=np.float32)  # expects 3
    with pytest.raises(ValueError, match="sh must be empty when sh_degree == 0"):
        cloud.sh_degree = 0
        cloud.sh = np.zeros(3, dtype=np.float32)
    cloud.sh_degree = 2
    # First failure should be the multiple check; then check exact length when num_points>0
    with pytest.raises(ValueError, match="sh length must be a multiple of 24, got 45"):
        cloud.sh = np.zeros(45, dtype=np.float32)  # expects multiple of 24 for degree 2
    with pytest.raises(ValueError, match=r"sh length must equal num_points \* \(\(sh_degree\+1\)\^2 - 1\) \* 3"):
        cloud.sh = np.zeros(48, dtype=np.float32)  # multiple of 24 but wrong for 1 point


def test_error_handling_file_operations():
    """Test error handling for file operations."""
    cloud = make_test_gaussian_cloud(include_sh=False)
    
    # Test saving to invalid path
    invalid_path = "/invalid/path/that/does/not/exist/test.spz"
    result = spz.save_spz(cloud, spz.PackOptions(), invalid_path)
    assert result is False
    
    # Test loading non-existent file - returns empty cloud instead of raising exception
    loaded_cloud = spz.load_spz("non_existent_file.spz", spz.UnpackOptions())
    assert loaded_cloud.num_points == 0
    assert loaded_cloud.sh_degree == 0
    
    # Test loading invalid file format - should return empty cloud
    invalid_file = os.path.join(tempfile.gettempdir(), "invalid.spz")
    with open(invalid_file, 'w') as f:
        f.write("This is not a valid SPZ file")
    
    loaded_cloud = spz.load_spz(invalid_file, spz.UnpackOptions())
    assert loaded_cloud.num_points == 0
    assert loaded_cloud.sh_degree == 0


def test_io_consistency_spz_format():
    """Test that SPZ format maintains consistency across save/load cycles."""
    original = make_test_gaussian_cloud(include_sh=True)
    filename = os.path.join(tempfile.gettempdir(), "consistency_test.spz")
    
    # Save and load multiple times
    for i in range(3):
        assert spz.save_spz(original, spz.PackOptions(), filename) is True
        loaded = spz.load_spz(filename, spz.UnpackOptions())
        
        # Basic properties should be identical
        assert loaded.num_points == original.num_points
        assert loaded.sh_degree == original.sh_degree
        assert loaded.antialiased == original.antialiased
        
        # Arrays should be close (accounting for compression)
        np.testing.assert_allclose(loaded.positions, original.positions, atol=1/2048.0)
        np.testing.assert_allclose(loaded.alphas, original.alphas, atol=0.01)
        
        # Use the loaded cloud as input for the next iteration
        original = loaded


def test_io_consistency_ply_format():
    """Test that PLY format maintains consistency and compatibility."""
    original = make_test_gaussian_cloud(include_sh=True)
    filename = os.path.join(tempfile.gettempdir(), "consistency_test.ply")
    
    # Save and load
    assert spz.save_splat_to_ply(original, spz.PackOptions(), filename) is True
    loaded = spz.load_splat_from_ply(filename, spz.UnpackOptions())
    
    # PLY format should preserve exact values (no compression)
    assert loaded.num_points == original.num_points
    assert loaded.sh_degree == original.sh_degree
    np.testing.assert_array_equal(loaded.positions, original.positions)
    np.testing.assert_array_equal(loaded.scales, original.scales)
    np.testing.assert_array_equal(loaded.rotations, original.rotations)
    np.testing.assert_array_equal(loaded.alphas, original.alphas)
    np.testing.assert_array_equal(loaded.colors, original.colors)
    np.testing.assert_array_equal(loaded.sh, original.sh)