"""Tests for spherical harmonics encoding and degrees."""

import os
import tempfile

import numpy as np
import pytest

import spz
from test_utils import sh_epsilon


def test_sh_encoding_for_zeros_and_edges():
    """Test spherical harmonics encoding for edge values and zeros."""
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


def test_spherical_harmonics_degree_coefficients():
    """Test that spherical harmonics coefficients match expected counts for different degrees."""
    cloud = spz.GaussianCloud()

    # Test degree 0 (no SH coefficients)
    cloud.sh_degree = 0
    cloud.sh = np.array([], dtype=np.float32)
    assert len(cloud.sh) == 0

    # Test degree 1 (9 coefficients per point: 3 coeffs × 3 channels)
    cloud.sh_degree = 1
    cloud.sh = np.array([0.0] * 9, dtype=np.float32)
    assert len(cloud.sh) == 9

    # Test degree 2 (24 coefficients per point: 8 coeffs × 3 channels)
    cloud.sh_degree = 2
    cloud.sh = np.array([0.0] * 24, dtype=np.float32)
    assert len(cloud.sh) == 24

    # Test degree 3 (45 coefficients per point: 15 coeffs × 3 channels)
    cloud.sh_degree = 3
    cloud.sh = np.array([0.0] * 45, dtype=np.float32)
    assert len(cloud.sh) == 45

    # Test multiple points
    cloud.positions = np.zeros(2 * 3, dtype=np.float32)
    cloud.sh_degree = 1
    cloud.sh = np.array([0.0] * 18, dtype=np.float32)  # 2 points × 9 coeffs
    assert len(cloud.sh) == 18


def test_spherical_harmonics_coordinate_transformation():
    """Test that spherical harmonics are properly transformed during coordinate conversion."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1
    cloud.antialiased = False

    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)

    original_sh = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], dtype=np.float32)
    cloud.sh = original_sh.copy()

    filename = os.path.join(tempfile.gettempdir(), "sh_coord_test.spz")

    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RUB
    assert spz.save_spz(cloud, pack_opts, filename) is True

    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RDF
    loaded_cloud = spz.load_spz(filename, unpack_opts)

    assert loaded_cloud.num_points == 1
    assert loaded_cloud.sh_degree == 1
    assert len(loaded_cloud.sh) == 9

    # SH coefficients should be transformed (not equal to original)
    assert not np.allclose(loaded_cloud.sh, original_sh, rtol=0, atol=1e-6)

    # Verify the transformation is consistent
    filename2 = os.path.join(tempfile.gettempdir(), "sh_coord_test2.spz")
    pack_opts2 = spz.PackOptions()
    pack_opts2.from_coord = spz.RDF
    assert spz.save_spz(loaded_cloud, pack_opts2, filename2) is True

    unpack_opts2 = spz.UnpackOptions()
    unpack_opts2.to_coord = spz.RDF
    loaded_cloud2 = spz.load_spz(filename2, unpack_opts2)

    np.testing.assert_array_almost_equal(loaded_cloud.sh, loaded_cloud2.sh, decimal=4)


# -----------------------------------------------------------------------------
# SH Degree 4 Tests
# -----------------------------------------------------------------------------


def test_sh_degree_4_coefficients():
    """Test that SH degree 4 produces correct number of coefficients."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 4

    # Degree 4: ((degree+1)^2 - 1) * 3 = (25 - 1) * 3 = 72 per point
    expected_coeffs_per_point = 72

    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.sh = np.zeros(expected_coeffs_per_point, dtype=np.float32)

    assert len(cloud.sh) == expected_coeffs_per_point
    assert cloud.sh_degree == 4


def test_sh_degree_4_save_load():
    """Test saving and loading SPZ files with SH degree 4."""
    num_points = 5
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 4

    rng = np.random.default_rng(42)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 72).astype(np.float32)

    filename = os.path.join(tempfile.gettempdir(), "sh_degree_4_test.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), filename) is True

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.num_points == num_points
    assert loaded.sh_degree == 4
    assert len(loaded.sh) == num_points * 72

    np.testing.assert_allclose(loaded.positions, cloud.positions, atol=1/2048.0)
    np.testing.assert_allclose(loaded.scales, cloud.scales, atol=1/16.0)
    np.testing.assert_allclose(loaded.alphas, cloud.alphas, atol=0.01)
    np.testing.assert_allclose(loaded.sh, cloud.sh, atol=sh_epsilon(4))


def test_sh_degree_4_ply_roundtrip():
    """Test PLY format with SH degree 4 preserves exact values."""
    num_points = 3
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 4

    rng = np.random.default_rng(123)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 72).astype(np.float32)

    filename = os.path.join(tempfile.gettempdir(), "sh_degree_4_test.ply")
    assert spz.save_splat_to_ply(cloud, spz.PackOptions(), filename) is True

    loaded = spz.load_splat_from_ply(filename, spz.UnpackOptions())
    assert loaded.num_points == num_points
    assert loaded.sh_degree == 4
    assert len(loaded.sh) == num_points * 72

    # PLY preserves exact float32 values
    np.testing.assert_allclose(loaded.positions, cloud.positions, rtol=0, atol=0)
    np.testing.assert_allclose(loaded.scales, cloud.scales, rtol=0, atol=0)
    np.testing.assert_allclose(loaded.sh, cloud.sh, rtol=0, atol=0)


def test_sh_degree_4_invalid_length():
    """Test that incorrect SH array length for degree 4 raises an error."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 4
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)

    with pytest.raises(ValueError, match="sh length must be a multiple of 72"):
        cloud.sh = np.zeros(45, dtype=np.float32)

    with pytest.raises(ValueError, match="sh length must equal num_points"):
        cloud.sh = np.zeros(144, dtype=np.float32)

