"""Tests for error handling and edge cases."""

import os
import tempfile

import numpy as np
import pytest

import spz
from test_utils import make_test_gaussian_cloud


def test_edge_cases_empty_arrays():
    """Test handling of empty arrays and edge cases."""
    cloud = spz.GaussianCloud()

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


def test_edge_cases_empty_cloud():
    """Test operations on empty GaussianCloud."""
    cloud = spz.GaussianCloud()

    # Test median_volume on empty cloud
    np.testing.assert_almost_equal(cloud.median_volume(), 0.01, decimal=5)

    # Test rotate_180_deg_about_x on empty cloud
    cloud.rotate_180_deg_about_x()

    # Test saving empty cloud
    filename = os.path.join(tempfile.gettempdir(), "empty_cloud.spz")
    result = spz.save_spz(cloud, spz.PackOptions(), filename)
    assert result is True

    # Test loading empty cloud
    loaded_cloud = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded_cloud.num_points == 0
    assert loaded_cloud.sh_degree == 0
    assert len(loaded_cloud.positions) == 0


def test_error_handling_invalid_arrays():
    """Test error handling for invalid array inputs."""
    cloud = spz.GaussianCloud()

    # Test with complex arrays that can't be converted to float32
    with pytest.raises(TypeError):
        complex_array = np.array([1+2j, 3+4j], dtype=np.complex64)
        cloud.positions = complex_array

    # Shape validations for properties
    cloud.sh_degree = 2
    cloud.positions = np.zeros(3, dtype=np.float32)  # one point
    with pytest.raises(ValueError, match="positions length must be a multiple of 3"):
        cloud.positions = np.array([1, 2], dtype=np.float32)
    with pytest.raises(ValueError, match=r"scales length must equal num_points \* 3"):
        cloud.scales = np.zeros(6, dtype=np.float32)
    with pytest.raises(ValueError, match=r"rotations length must equal num_points \* 4"):
        cloud.rotations = np.zeros(8, dtype=np.float32)
    with pytest.raises(ValueError, match=r"colors length must equal num_points \* 3"):
        cloud.colors = np.zeros(6, dtype=np.float32)
    with pytest.raises(ValueError, match="sh must be empty when sh_degree == 0"):
        cloud.sh_degree = 0
        cloud.sh = np.zeros(3, dtype=np.float32)
    cloud.sh_degree = 2
    with pytest.raises(ValueError, match="sh length must be a multiple of 24, got 45"):
        cloud.sh = np.zeros(45, dtype=np.float32)
    with pytest.raises(ValueError, match=r"sh length must equal num_points \* \(\(sh_degree\+1\)\^2 - 1\) \* 3"):
        cloud.sh = np.zeros(48, dtype=np.float32)


def test_error_handling_file_operations():
    """Test error handling for file operations."""
    cloud = make_test_gaussian_cloud(include_sh=False)

    # Test saving to invalid path
    invalid_path = "/invalid/path/that/does/not/exist/test.spz"
    result = spz.save_spz(cloud, spz.PackOptions(), invalid_path)
    assert result is False

    # Test loading non-existent file - returns empty cloud
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

