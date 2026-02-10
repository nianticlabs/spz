"""Tests for SH quantization bit rates."""

import os
import tempfile

import numpy as np
import pytest

import spz

# Skip all tests in this module if extension support is not available
pytestmark = pytest.mark.skipif(
    not spz.has_extension_support(),
    reason="SPZ_BUILD_EXTENSIONS is not enabled"
)


def test_pack_options_sh_bits_properties():
    """Test that PackOptions has sh1_bits and sh_rest_bits properties."""
    opts = spz.PackOptions()

    # Test default values (from C++ constants: DEFAULT_SH1_BITS=5, DEFAULT_SH_REST_BITS=4)
    assert opts.sh1_bits == 5
    assert opts.sh_rest_bits == 4

    # Test setting custom values
    opts.sh1_bits = 8
    opts.sh_rest_bits = 6

    assert opts.sh1_bits == 8
    assert opts.sh_rest_bits == 6


def test_sh_quantization_8bit():
    """Test SH quantization with maximum 8-bit precision."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 2
    cloud.positions = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    cloud.scales = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.alphas = np.array([0.0], dtype=np.float32)
    cloud.colors = np.array([0.0, 0.0, 0.0], dtype=np.float32)

    original_sh = np.array([
        0.1, 0.2, 0.3,  # coeff 0 (degree 1)
        0.4, 0.5, 0.6,  # coeff 1 (degree 1)
        0.7, 0.8, 0.9,  # coeff 2 (degree 1)
        -0.1, -0.2, -0.3,  # coeff 3 (degree 2)
        -0.4, -0.5, -0.6,  # coeff 4 (degree 2)
        0.15, 0.25, 0.35,  # coeff 5 (degree 2)
        0.45, 0.55, 0.65,  # coeff 6 (degree 2)
        0.75, 0.85, 0.95,  # coeff 7 (degree 2)
    ], dtype=np.float32)
    cloud.sh = original_sh.copy()

    opts = spz.PackOptions()
    opts.sh1_bits = 8
    opts.sh_rest_bits = 8

    filename = os.path.join(tempfile.gettempdir(), "sh_8bit_test.spz")
    assert spz.save_spz(cloud, opts, filename) is True

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    np.testing.assert_allclose(loaded.sh, original_sh, atol=0.01)


def test_sh_quantization_low_bits():
    """Test SH quantization with lower bit precision."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1
    cloud.positions = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    cloud.scales = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.alphas = np.array([0.0], dtype=np.float32)
    cloud.colors = np.array([0.0, 0.0, 0.0], dtype=np.float32)

    original_sh = np.array([0.5, 0.5, 0.5, -0.5, -0.5, -0.5, 0.25, 0.25, 0.25], dtype=np.float32)
    cloud.sh = original_sh.copy()

    opts = spz.PackOptions()
    opts.sh1_bits = 1
    opts.sh_rest_bits = 1

    filename = os.path.join(tempfile.gettempdir(), "sh_lowbit_test.spz")
    assert spz.save_spz(cloud, opts, filename) is True

    loaded = spz.load_spz(filename, spz.UnpackOptions())

    assert len(loaded.sh) == 9
    # Large tolerance expected with 1-bit quantization
    np.testing.assert_allclose(loaded.sh, original_sh, atol=1.0)


def test_sh_bits_comparison():
    """Test that higher bit rates produce better precision."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 2

    num_points = 10
    rng = np.random.default_rng(42)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 24).astype(np.float32)

    errors = {}
    for bits in [4, 6, 8]:
        opts = spz.PackOptions()
        opts.sh1_bits = bits
        opts.sh_rest_bits = bits

        filename = os.path.join(tempfile.gettempdir(), f"sh_bits_{bits}.spz")
        assert spz.save_spz(cloud, opts, filename) is True

        loaded = spz.load_spz(filename, spz.UnpackOptions())
        errors[bits] = np.abs(loaded.sh - cloud.sh).max()

    # Higher bits should give lower error
    assert errors[8] <= errors[6] <= errors[4]


def test_sh_degree_4_with_custom_bits():
    """Test SH degree 4 with custom bit rates."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 4

    num_points = 5
    rng = np.random.default_rng(42)
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(0.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-1.0, 1.0, size=num_points * 72).astype(np.float32)

    opts = spz.PackOptions()
    opts.version = 3
    opts.sh1_bits = 6
    opts.sh_rest_bits = 3

    filename = os.path.join(tempfile.gettempdir(), "sh4_custom_bits.spz")
    assert spz.save_spz(cloud, opts, filename) is True

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.sh_degree == 4
    assert len(loaded.sh) == num_points * 72

    # Degree 1 coefficients should have reasonable error with 6-bit quantization
    np.testing.assert_allclose(loaded.sh[:9], cloud.sh[:9], rtol=0, atol=0.02)

