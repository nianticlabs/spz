"""Tests for SPZ version compatibility and regression."""

import os
import tempfile

import numpy as np

import spz
from test_utils import times


def test_pack_options_version_property():
    """Test that PackOptions has version property."""
    opts = spz.PackOptions()

    # Default should be latest version (3)
    assert opts.version == 3

    opts.version = 2
    assert opts.version == 2

    opts.version = 1
    assert opts.version == 1


def test_spz_version_3_quaternion_encoding():
    """Test SPZ version 3 with smallest-three quaternion encoding."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1

    num_points = 5
    rng = np.random.default_rng(42)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 9).astype(np.float32)

    opts = spz.PackOptions()
    opts.version = 3

    filename = os.path.join(tempfile.gettempdir(), "version_3_test.spz")
    assert spz.save_spz(cloud, opts, filename) is True

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.num_points == num_points

    # Verify quaternions are normalized and preserve rotation
    for i in range(num_points):
        q_orig = cloud.rotations[i*4:(i+1)*4]
        q_loaded = loaded.rotations[i*4:(i+1)*4]

        assert abs(np.linalg.norm(q_loaded) - 1.0) < 1e-4

        q_orig_norm = q_orig / np.linalg.norm(q_orig)
        test_vec = np.array([1.0, 2.0, 3.0])

        rotated_orig = times(q_orig_norm, test_vec)
        rotated_loaded = times(q_loaded, test_vec)

        cosine = np.dot(rotated_orig, rotated_loaded) / (np.linalg.norm(rotated_orig) * np.linalg.norm(rotated_loaded))
        assert abs(cosine - 1.0) < 1e-3


def test_spz_version_2_first_three_encoding():
    """Test SPZ version 2 with first-three quaternion encoding (legacy)."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1

    num_points = 5
    rng = np.random.default_rng(42)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 9).astype(np.float32)

    opts = spz.PackOptions()
    opts.version = 2

    filename = os.path.join(tempfile.gettempdir(), "version_2_test.spz")
    assert spz.save_spz(cloud, opts, filename) is True

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.num_points == num_points

    np.testing.assert_allclose(loaded.positions, cloud.positions, atol=1/2048.0)
    np.testing.assert_allclose(loaded.alphas, cloud.alphas, atol=0.01)


def test_spz_version_compatibility_positions():
    """Test that position encoding is consistent across versions 2 and 3."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0

    cloud.positions = np.array([
        1.5, -2.5, 3.5,
        0.001, 0.002, 0.003,
        100.0, -100.0, 50.0
    ], dtype=np.float32)
    cloud.scales = np.zeros(9, dtype=np.float32)
    cloud.rotations = np.array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1], dtype=np.float32)
    cloud.alphas = np.array([0.5, 0.5, 0.5], dtype=np.float32)
    cloud.colors = np.zeros(9, dtype=np.float32)
    cloud.sh = np.array([], dtype=np.float32)

    for version in [2, 3]:
        opts = spz.PackOptions()
        opts.version = version

        filename = os.path.join(tempfile.gettempdir(), f"version_{version}_positions.spz")
        assert spz.save_spz(cloud, opts, filename) is True

        loaded = spz.load_spz(filename, spz.UnpackOptions())

        np.testing.assert_allclose(loaded.positions, cloud.positions, atol=1/2048.0,
                                   err_msg=f"Position mismatch for version {version}")


def test_spz_version_3_sh_degree_4():
    """Test that version 3 correctly handles SH degree 4."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 4

    num_points = 3
    rng = np.random.default_rng(42)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 72).astype(np.float32)

    opts = spz.PackOptions()
    opts.version = 3

    filename = os.path.join(tempfile.gettempdir(), "version_3_sh4.spz")
    assert spz.save_spz(cloud, opts, filename) is True

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.sh_degree == 4
    assert len(loaded.sh) == num_points * 72


def test_version_roundtrip_consistency():
    """Test that saving with a version and loading produces consistent results."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 2
    cloud.antialiased = True

    num_points = 10
    rng = np.random.default_rng(123)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 24).astype(np.float32)

    for version in [2, 3]:
        opts = spz.PackOptions()
        opts.version = version

        current = cloud
        for iteration in range(3):
            filename = os.path.join(tempfile.gettempdir(), f"roundtrip_v{version}_iter{iteration}.spz")
            assert spz.save_spz(current, opts, filename) is True
            loaded = spz.load_spz(filename, spz.UnpackOptions())

            assert loaded.num_points == cloud.num_points
            assert loaded.sh_degree == cloud.sh_degree
            assert loaded.antialiased == cloud.antialiased

            current = loaded

        np.testing.assert_allclose(current.positions, cloud.positions, atol=1/1024.0,
                                   err_msg=f"Position drift after roundtrips for version {version}")


def test_version_3_with_all_features():
    """Test version 3 with SH degree 4, custom bits."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 4
    cloud.antialiased = True

    num_points = 20
    rng = np.random.default_rng(999)
    cloud.positions = rng.uniform(-10.0, 10.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-3.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(-2.0, 2.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(-0.5, 1.5, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-1.5, 1.5, size=num_points * 72).astype(np.float32)

    opts = spz.PackOptions()
    opts.version = 3
    opts.sh1_bits = 8
    opts.sh_rest_bits = 6
    opts.from_coord = spz.RUB

    filename = os.path.join(tempfile.gettempdir(), "all_features.spz")
    assert spz.save_spz(cloud, opts, filename) is True

    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RUB

    loaded = spz.load_spz(filename, unpack_opts)

    assert loaded.num_points == num_points
    assert loaded.sh_degree == 4
    assert loaded.antialiased
    assert len(loaded.sh) == num_points * 72

    np.testing.assert_allclose(loaded.sh, cloud.sh, rtol=0, atol=0.05)

