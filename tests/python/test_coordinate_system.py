"""Tests for coordinate system conversions and transformations."""

import os
import tempfile

import numpy as np

import spz


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
    assert hasattr(spz, 'LBD')
    assert hasattr(spz, 'RBD')
    assert hasattr(spz, 'LBU')
    assert hasattr(spz, 'RBU')
    assert hasattr(spz, 'LFD')
    assert hasattr(spz, 'RFD')
    assert hasattr(spz, 'LFU')
    assert hasattr(spz, 'RFU')

    # Test that all enum values are unique
    enum_values = [
        spz.UNSPECIFIED, spz.LDB, spz.RDB, spz.LUB, spz.RUB,
        spz.LDF, spz.RDF, spz.LUF, spz.RUF,
        spz.LBD, spz.RBD, spz.LBU, spz.RBU,
        spz.LFD, spz.RFD, spz.LFU, spz.RFU,
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
    expected_pos = np.array([1.0, -2.0, -3.0], dtype=np.float32)
    np.testing.assert_allclose(loaded_cloud.positions, expected_pos, atol=1/2048.0)

    # Quaternion Y and Z components should be flipped (but not W)
    expected_rot = np.array([0.1, -0.2, -0.3, 0.9], dtype=np.float32)
    loaded_rot_norm = loaded_cloud.rotations / np.linalg.norm(loaded_cloud.rotations)
    expected_rot_norm = expected_rot / np.linalg.norm(expected_rot)
    np.testing.assert_allclose(loaded_rot_norm, expected_rot_norm, atol=1e-3)

    # Test conversion from RDF to LUF (more complex transformation)
    pack_opts.from_coord = spz.RDF
    unpack_opts.to_coord = spz.LUF

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

    # RDF to LUF should flip X and Y
    expected_pos2 = np.array([-1.0, -2.0, 3.0], dtype=np.float32)
    np.testing.assert_allclose(loaded_cloud2.positions, expected_pos2, atol=1/2048.0)


def test_convert_coordinates_method():
    """Directly test in-place coordinate conversion binding."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.rotations = np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32)

    # Convert to RDF: should flip Y and Z for positions and quaternion components
    cloud.convert_coordinates(spz.RUB, spz.RDF)
    np.testing.assert_allclose(cloud.positions, np.array([1.0, -2.0, -3.0], dtype=np.float32), atol=1e-6)
    np.testing.assert_allclose(cloud.rotations, np.array([0.1, -0.2, -0.3, 0.9], dtype=np.float32), atol=1e-6)

    # Convert back to RUB: should restore original values
    cloud.convert_coordinates(spz.RDF, spz.RUB)
    np.testing.assert_allclose(cloud.positions, np.array([1.0, 2.0, 3.0], dtype=np.float32), atol=1e-6)
    np.testing.assert_allclose(cloud.rotations, np.array([0.1, 0.2, 0.3, 0.9], dtype=np.float32), atol=1e-6)

