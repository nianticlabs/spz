"""Tests for GaussianCloud class properties and methods."""

import numpy as np
import pytest

import spz


def test_gaussian_cloud_initialization():
    """Test GaussianCloud initialization and basic properties."""
    cloud = spz.GaussianCloud()

    # Test default values
    assert cloud.num_points == 0
    assert cloud.sh_degree == 0
    assert not cloud.antialiased

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
    assert cloud.antialiased

    # Test setting array properties
    positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.positions = positions
    np.testing.assert_allclose(cloud.positions, positions, rtol=0, atol=0)
    assert cloud.num_points == 1

    scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.scales = scales
    np.testing.assert_allclose(cloud.scales, scales, rtol=0, atol=0)

    rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.rotations = rotations
    np.testing.assert_allclose(cloud.rotations, rotations, rtol=0, atol=0)

    alphas = np.array([0.5], dtype=np.float32)
    cloud.alphas = alphas
    np.testing.assert_allclose(cloud.alphas, alphas, rtol=0, atol=0)

    colors = np.array([1.0, 0.0, 0.0], dtype=np.float32)
    cloud.colors = colors
    np.testing.assert_allclose(cloud.colors, colors, rtol=0, atol=0)

    # For sh_degree=2, per-point SH length must be 24
    sh = np.zeros(24, dtype=np.float32)
    cloud.sh = sh
    np.testing.assert_allclose(cloud.sh, sh, rtol=0, atol=0)


def test_gaussian_cloud_array_dtype_handling():
    """Test that GaussianCloud properly handles different array dtypes."""
    cloud = spz.GaussianCloud()

    # Test with different numpy dtypes - should convert to float32
    for dtype in [np.float64, np.int32, np.float32]:
        data = np.array([1.0, 2.0, 3.0], dtype=dtype)
        cloud.positions = data
        # The returned array should always be float32
        assert cloud.positions.dtype == np.float32
        np.testing.assert_allclose(cloud.positions, [1.0, 2.0, 3.0], rtol=0, atol=0)

    # Test that non-numeric types are rejected
    with pytest.raises(TypeError, match="incompatible function arguments"):
        cloud.positions = np.array(['a', 'b', 'c'], dtype=np.str_)

    # Test that complex types are rejected
    with pytest.raises(TypeError, match="incompatible function arguments"):
        cloud.positions = np.array([1+2j, 3+4j, 5+6j], dtype=np.complex64)


def test_gaussian_cloud_rotate_180_deg_about_x():
    """Test the rotate_180_deg_about_x method."""
    cloud = spz.GaussianCloud()

    original_pos = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    original_rot = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)

    cloud.positions = original_pos.copy()
    cloud.rotations = original_rot.copy()

    # Apply rotation - converts from RUB to RDF coordinates
    cloud.rotate_180_deg_about_x()

    # Check that positions have been transformed (Y and Z should be flipped)
    expected_pos = np.array([1.0, -2.0, -3.0], dtype=np.float32)
    np.testing.assert_array_almost_equal(cloud.positions, expected_pos, decimal=5)

    # Check that rotations have been transformed
    expected_rot = np.array([0.1, -0.2, -0.3, 0.9], dtype=np.float32)
    np.testing.assert_array_almost_equal(cloud.rotations, expected_rot, decimal=5)


def test_gaussian_cloud_median_volume():
    """Test the median_volume method calculates volume correctly."""
    cloud = spz.GaussianCloud()

    # Test empty cloud - should return 0.01 according to C++ implementation
    np.testing.assert_almost_equal(cloud.median_volume(), 0.01, decimal=5)

    # Test with actual points (3 points)
    cloud.positions = np.zeros(3 * 3, dtype=np.float32)
    cloud.scales = np.array([
        -1.0, -1.0, -1.0,  # First gaussian: scale sum = -3
        0.0, 0.0, 0.0,     # Second gaussian: scale sum = 0
        1.0, 1.0, 1.0      # Third gaussian: scale sum = 3
    ], dtype=np.float32)

    median_vol = cloud.median_volume()
    expected_vol = (4.0 / 3.0) * np.pi * np.exp(0.0)
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
    expected_vol = (4.0 / 3.0) * np.pi * np.exp(0.0)
    np.testing.assert_almost_equal(median_vol, expected_vol, decimal=5)


def test_quaternion_normalization_during_packing():
    """Test that quaternions are normalized during SPZ packing/unpacking."""
    import os
    import tempfile

    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.antialiased = False

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

    filename = os.path.join(tempfile.gettempdir(), "quat_normalization_test.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), filename) is True

    loaded_cloud = spz.load_spz(filename, spz.UnpackOptions())

    # Verify that the quaternions are now normalized
    assert len(loaded_cloud.rotations) == 8

    q1 = loaded_cloud.rotations[0:4]
    q1_norm = np.linalg.norm(q1)
    assert abs(q1_norm - 1.0) < 1e-4

    q2 = loaded_cloud.rotations[4:8]
    q2_norm = np.linalg.norm(q2)
    assert abs(q2_norm - 1.0) < 1e-4

    # Verify the orientation is preserved
    expected_q1 = non_normalized_quats[0:4] / np.linalg.norm(non_normalized_quats[0:4])
    expected_q2 = non_normalized_quats[4:8] / np.linalg.norm(non_normalized_quats[4:8])

    def quaternions_equivalent(q1, q2, tolerance=1e-2):
        return (np.allclose(q1, q2, atol=tolerance) or
                np.allclose(q1, -q2, atol=tolerance))

    assert quaternions_equivalent(q1, expected_q1) or quaternions_equivalent(q1, -expected_q1)
    assert quaternions_equivalent(q2, expected_q2) or quaternions_equivalent(q2, -expected_q2)

