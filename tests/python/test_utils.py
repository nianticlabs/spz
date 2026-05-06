"""Shared test utilities, constants, and helper functions."""

import numpy as np
from scipy.spatial.transform import Rotation

import spz

# -----------------------------------------------------------------------------
# Epsilon constants (from the C++ tests)
# -----------------------------------------------------------------------------


def sh_epsilon(bits):
    """
    Calculate the maximum quantization error (epsilon) for SH quantization with given bit precision.

    For N-bit quantization:
    - Number of buckets = 2^N
    - Half bucket size (max quantization error) = 1.0 / (2^N) = 1.0 / (1 << bits)
    - Rounding error from quantization step = 0.5 / 128.0 (constant)
    - Total max error = half_bucket + rounding_error

    Args:
        bits: Number of quantization bits (1-8)

    Returns:
        Maximum quantization error (epsilon) for the given bit precision
    """
    return 1.0 / (1 << bits) + 0.5 / 128.0

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
    SciPy's Rotation.from_rotvec returns a quaternion in [x, y, z, w] order.
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


def make_single_point_cloud():
    """Return a 1-point GaussianCloud with simple data at (1, 2, 3)."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.antialiased = False
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.sh = np.array([], dtype=np.float32)
    return cloud


def sh_plus_rotation_matrix(band_idx):
    """Return the (2l+1)×(2l+1) SH rotation matrix for R_x(+π/2) on band l=band_idx+1."""
    l = band_idx + 1
    n = 2 * l + 1
    coeffs_before = band_idx * (band_idx + 2)  # == C++ bandStart = band*(band+2)
    total_coeffs = coeffs_before + n
    sh_size = total_coeffs * 3

    matrix = np.zeros((n, n))
    for col in range(n):
        sh = np.zeros(sh_size, dtype=np.float32)
        sh[(coeffs_before + col) * 3] = 1.0  # unit vector in R channel of coefficient col
        cloud = spz.GaussianCloud()
        cloud.sh_degree = l
        cloud.positions = np.array([1.0, 0.0, 0.0], dtype=np.float32)
        cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
        cloud.sh = sh
        cloud.convert_coordinates(spz.RUF, spz.RBU)
        for row in range(n):
            matrix[row, col] = cloud.sh[(coeffs_before + row) * 3]
    return matrix


def make_test_gaussian_cloud(include_sh):
    """
    Create a GaussianCloud with two splats for testing.
    If include_sh is True then spherical harmonics (SH) are added and sh_degree is set to 3.
    """
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

