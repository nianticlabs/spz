"""Tests for SPZ and PLY file I/O operations."""

import os
import tempfile

import numpy as np
import pytest

import spz
from test_utils import (
    sh_epsilon,
    make_test_gaussian_cloud,
    normalized,
    read_file,
    times,
)


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
    np.testing.assert_allclose(dst.sh, src.sh, atol=sh_epsilon(4))
    # Check degree‑1 SH (first 9 coefficients) with extra precision.
    np.testing.assert_allclose(dst.sh[0:9], src.sh[0:9], atol=sh_epsilon(5))
    np.testing.assert_allclose(dst.sh[45:45 + 9], src.sh[45:45 + 9], atol=sh_epsilon(5))


def test_save_load_packed_format_large_splat():
    """Test saving and loading large SPZ files with many points."""
    num_points = 50000
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
    # PLY format preserves exact float32 values
    np.testing.assert_allclose(dst.positions, src.positions, rtol=0, atol=0)
    np.testing.assert_allclose(dst.scales, src.scales, rtol=0, atol=0)
    np.testing.assert_allclose(dst.rotations, src.rotations, rtol=0, atol=0)
    np.testing.assert_allclose(dst.alphas, src.alphas, rtol=0, atol=0)
    np.testing.assert_allclose(dst.colors, src.colors, rtol=0, atol=0)
    if include_sh:
        np.testing.assert_allclose(dst.sh, src.sh, rtol=0, atol=0)
    else:
        assert len(dst.sh) == 0


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
    np.testing.assert_allclose(loaded.positions, original.positions, rtol=0, atol=0)
    np.testing.assert_allclose(loaded.scales, original.scales, rtol=0, atol=0)
    np.testing.assert_allclose(loaded.rotations, original.rotations, rtol=0, atol=0)
    np.testing.assert_allclose(loaded.alphas, original.alphas, rtol=0, atol=0)
    np.testing.assert_allclose(loaded.colors, original.colors, rtol=0, atol=0)
    np.testing.assert_allclose(loaded.sh, original.sh, rtol=0, atol=0)


def test_compression_precision_validation():
    """Test that compression maintains expected precision levels."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1
    cloud.antialiased = False

    # Test with specific values that test the compression boundaries
    cloud.positions = np.array([1.0, -1.0, 0.5], dtype=np.float32)
    cloud.scales = np.array([1.0, -1.0, 0.5], dtype=np.float32)
    cloud.rotations = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.5, -0.5, 0.25], dtype=np.float32)
    cloud.sh = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], dtype=np.float32)

    filename = os.path.join(tempfile.gettempdir(), "compression_precision_test.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), filename) is True

    loaded_cloud = spz.load_spz(filename, spz.UnpackOptions())

    # Verify precision according to C++ implementation
    np.testing.assert_allclose(loaded_cloud.positions, cloud.positions, atol=1/2048.0)
    np.testing.assert_allclose(loaded_cloud.scales, cloud.scales, atol=1/32.0)

    # Quaternions should be normalized and preserve orientation
    loaded_quat = loaded_cloud.rotations
    assert abs(np.linalg.norm(loaded_quat) - 1.0) < 1e-4

    np.testing.assert_allclose(loaded_cloud.alphas, cloud.alphas, atol=0.01)
    np.testing.assert_allclose(loaded_cloud.colors, cloud.colors, atol=0.01)
    np.testing.assert_allclose(loaded_cloud.sh[0:9], cloud.sh[0:9], atol=sh_epsilon(5))

    # Verify the function matches the C++ implementation
    # Note: The rounding error uses 128.0 because quantization maps [-1,1] to [0,256]
    assert sh_epsilon(4) == 2.0 / 32.0 + 0.5 / 128.0
    assert sh_epsilon(5) == 2.0 / 64.0 + 0.5 / 128.0


def test_performance_large_cloud():
    """Test performance with a large cloud and verify timing."""
    import time

    num_points = 10000
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 2

    rng = np.random.default_rng(42)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 24).astype(np.float32)

    filename = os.path.join(tempfile.gettempdir(), "large_cloud_performance.spz")
    start_time = time.time()
    result = spz.save_spz(cloud, spz.PackOptions(), filename)
    save_time = time.time() - start_time

    assert result is True
    assert save_time < 5.0

    start_time = time.time()
    loaded_cloud = spz.load_spz(filename, spz.UnpackOptions())
    load_time = time.time() - start_time

    assert loaded_cloud.num_points == num_points
    assert load_time < 5.0

