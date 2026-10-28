"""Tests for SPZ extension round-trip through NGSP v4 and legacy gzip formats."""

import gzip
import io
import os
import struct
import tempfile

import numpy as np
import pytest

import spz


def _make_cloud(num_points=5, seed=0):
    rng = np.random.default_rng(seed)
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(-1.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 9).astype(np.float32)
    return cloud


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_safe_orbit_extension_round_trip():
    """SafeOrbit extension survives a save/load round-trip"""
    cloud = _make_cloud()

    ext = spz.SpzExtensionSafeOrbitCameraAdobe()
    ext.safe_orbit_elevation_min = -0.5
    ext.safe_orbit_elevation_max = 1.2
    ext.safe_orbit_radius_min = 0.3
    cloud.extensions = [ext]

    filename = os.path.join(tempfile.gettempdir(), "ext_round_trip.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), filename) is True

    # FlagHasExtensions (0x2) must be set in the header flags byte (offset 14).
    with open(filename, "rb") as f:
        header = f.read(32)
    flags = header[14]
    assert flags & 0x2, f"FlagHasExtensions not set in saved file (flags=0x{flags:02x})"

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.num_points == cloud.num_points
    assert len(loaded.extensions) == 1

    loaded_ext = loaded.extensions[0]
    assert isinstance(loaded_ext, spz.SpzExtensionSafeOrbitCameraAdobe)
    assert abs(loaded_ext.safe_orbit_elevation_min - ext.safe_orbit_elevation_min) < 1e-5
    assert abs(loaded_ext.safe_orbit_elevation_max - ext.safe_orbit_elevation_max) < 1e-5
    assert abs(loaded_ext.safe_orbit_radius_min - ext.safe_orbit_radius_min) < 1e-5


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_no_extension_flag_when_no_extensions():
    """FlagHasExtensions must NOT be set when the cloud has no extensions."""
    cloud = _make_cloud()
    cloud.extensions = []

    filename = os.path.join(tempfile.gettempdir(), "no_ext.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), filename) is True

    with open(filename, "rb") as f:
        header = f.read(32)
    flags = header[14]
    assert not (flags & 0x2), f"FlagHasExtensions unexpectedly set (flags=0x{flags:02x})"


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_toc_byte_offset_advances_with_extensions():
    """tocByteOffset must be > 32 when extensions are correctly saved, 32 when absent."""
    cloud = _make_cloud()

    # Without extensions
    filename_no_ext = os.path.join(tempfile.gettempdir(), "toc_no_ext.spz")
    cloud.extensions = []
    assert spz.save_spz(cloud, spz.PackOptions(), filename_no_ext) is True
    with open(filename_no_ext, "rb") as f:
        hdr_no_ext = f.read(32)
    tbo_no_ext = struct.unpack_from("<I", hdr_no_ext, 16)[0]
    assert tbo_no_ext == 32, f"Expected tocByteOffset=32 with no extensions, got {tbo_no_ext}"

    # With extension
    ext = spz.SpzExtensionSafeOrbitCameraAdobe()
    cloud.extensions = [ext]
    filename_ext = os.path.join(tempfile.gettempdir(), "toc_with_ext.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), filename_ext) is True
    with open(filename_ext, "rb") as f:
        hdr_ext = f.read(32)
    tbo_ext = struct.unpack_from("<I", hdr_ext, 16)[0]
    assert tbo_ext > 32, f"Expected tocByteOffset>32 with extension, got {tbo_ext}"


def _attribute(name, values, **kwargs):
    return spz.SpzExtensionCustomAttributeAdobe(name, np.asarray(values, np.float32).ravel(), **kwargs)


def _custom_attributes(cloud):
    return [e for e in cloud.extensions if isinstance(e, spz.SpzExtensionCustomAttributeAdobe)]


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_custom_attributes_round_trip():
    """Float32 attributes load exactly, quantized ones within half a level, next to other extensions."""
    n = 200
    rng = np.random.default_rng(1)
    exact = rng.normal(size=(n, 2))
    coarse = rng.uniform(-5.0, 5.0, size=(n, 3)).astype(np.float32)
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_radius_min = 0.3
    cloud = _make_cloud(n)
    cloud.extensions = [_attribute("exact", exact, dim=2), orbit, _attribute("coarse", coarse, dim=3, bits=6)]

    loaded = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))

    assert loaded.num_points == n
    assert any(isinstance(e, spz.SpzExtensionSafeOrbitCameraAdobe) for e in loaded.extensions)
    attrs = _custom_attributes(loaded)
    assert [(a.name, a.dim, a.bits) for a in attrs] == [("exact", 2, 32), ("coarse", 3, 6)]
    np.testing.assert_array_equal(attrs[0].values, exact.astype(np.float32).ravel())
    half_level = (coarse.max(axis=0) - coarse.min(axis=0)) / 63 / 2
    assert np.all(np.abs(attrs[1].values.reshape(n, 3) - coarse) <= half_level * 1.001)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("from_coord,to_coord", [
    (spz.RDF, spz.RUB),
    (spz.RUB, spz.RBD),
    (spz.LUF, spz.RDF),
])
def test_vector3_attribute_converts_like_positions(from_coord, to_coord):
    """A Vector3 attribute holding the positions still equals them after any save/load conversion."""
    n = 50
    cloud = _make_cloud(n)
    positions = cloud.positions
    vectors = _attribute("p", positions, dim=3, kind=spz.SpzCustomAttributeKind.VECTOR3)
    scalars = _attribute("s", positions, dim=3)
    cloud.extensions = [vectors, scalars]
    pack_options = spz.PackOptions()
    pack_options.from_coord = from_coord
    unpack_options = spz.UnpackOptions()
    unpack_options.to_coord = to_coord

    loaded = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, pack_options), unpack_options)

    loaded_vectors, loaded_scalars = _custom_attributes(loaded)
    np.testing.assert_allclose(loaded_vectors.values, loaded.positions, atol=2.0 / 4096)
    np.testing.assert_array_equal(loaded_scalars.values, positions)
    np.testing.assert_array_equal(vectors.values, positions)  # saving didn't touch the caller's record


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("make_attributes", [
    lambda n: [_attribute("short", np.zeros(n - 1))],
    lambda n: [_attribute("", np.zeros(n))],
    lambda n: [_attribute("twice", np.zeros(n)), _attribute("twice", np.ones(n))],
    lambda n: [_attribute("bits", np.zeros(n), bits=9)],
    lambda n: [_attribute("vec", np.zeros(n * 2), dim=2, kind=spz.SpzCustomAttributeKind.VECTOR3)],
    lambda n: [_attribute("nan", np.full(n, np.nan), bits=8)],
], ids=["length", "empty-name", "duplicate-name", "bits", "vector-dim", "non-finite-quantized"])
def test_invalid_custom_attribute_fails_to_save(make_attributes):
    cloud = _make_cloud(5)
    cloud.extensions = make_attributes(5)
    with pytest.raises(RuntimeError):
        spz.save_spz_to_buffer(cloud, spz.PackOptions())


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("offset,value", [
    (43, struct.pack("<H", 2)),       # dim: data covers half the points it claims
    (43, struct.pack("<H", 65535)),   # dim: claims far more data than the frame can hold
    (46, bytes([9])),                 # bits: not a defined precision
], ids=["dim-mismatch", "dim-huge", "bits"])
def test_corrupt_custom_attribute_is_dropped(offset, value):
    """A malformed record is skipped; the core splat data still loads."""
    cloud = _make_cloud(5)
    cloud.extensions = [_attribute("a", np.arange(5))]
    data = bytearray(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    # header (32) + record type and length (8) + name length (2) + "a" (1) -> dim at 43, kind 45, bits 46
    data[offset:offset + len(value)] = value

    loaded = spz.load_spz_from_buffer(bytes(data))

    assert loaded.num_points == 5
    assert _custom_attributes(loaded) == []
