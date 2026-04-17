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
    """SafeOrbit extension survives a save/load round-trip through NGSP v4."""
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
    """tocByteOffset must be > 32 when extensions are present, 32 when absent."""
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


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_v3_gzip_extension_loading():
    """Extension data embedded in a legacy v3 gzip SPZ file is read correctly.

    Since the library no longer writes gzip files, we construct the bytes by hand
    using the documented legacy format so the read path can be tested independently
    of the write path.

    Legacy gzip layout (all little-endian, then gzip-compressed):
      LegacyPackedGaussiansHeader (16 bytes)
      positions  — numPoints × 9 bytes (24-bit fixed-point xyz)
      alphas     — numPoints × 1 byte
      colors     — numPoints × 3 bytes
      scales     — numPoints × 3 bytes
      rotations  — numPoints × 4 bytes (v3 = smallest-three quaternion)
      sh         — 0 bytes for sh_degree=0
      extensions — ILV records: type(u32) + byteLength(u32) + payload
    """
    NGSP_MAGIC = 0x5053474E
    FLAG_HAS_EXTENSIONS = 0x2
    SPZ_ADOBE_SAFE_ORBIT = 0xADBE0002

    elev_min, elev_max, radius_min = -0.5, 1.2, 0.3
    num_points = 1
    sh_degree = 0
    fractional_bits = 12

    # LegacyPackedGaussiansHeader: magic(4) version(4) numPoints(4) shDegree(1)
    #   fractionalBits(1) flags(1) reserved(1) — total 16 bytes
    header = struct.pack(
        "<IIIBBBx",
        NGSP_MAGIC, 3, num_points, sh_degree, fractional_bits, FLAG_HAS_EXTENSIONS,
    )

    # One point at the origin: all-zero 24-bit fixed-point positions
    positions = b"\x00" * (num_points * 9)
    alphas    = b"\x80" * num_points               # sigmoid^-1(128/255) ≈ 0
    colors    = b"\x80" * (num_points * 3)
    scales    = b"\x80" * (num_points * 3)

    # Identity quaternion in smallest-three encoding (v3+).
    # w=1 is the largest component (index 3); x=y=z=0.
    # Packed: iLargest=3 in top 2 bits, three zero 10-bit fields → 0xC0000000 LE.
    rotations = struct.pack("<I", 0xC0000000) * num_points

    # SafeOrbit extension: type(4) + byteLength(4) + 3 floats(12)
    payload   = struct.pack("<fff", elev_min, elev_max, radius_min)
    extension = struct.pack("<II", SPZ_ADOBE_SAFE_ORBIT, len(payload)) + payload

    raw = header + positions + alphas + colors + scales + rotations + extension

    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(raw)

    filename = os.path.join(tempfile.gettempdir(), "v3_ext_legacy.spz")
    with open(filename, "wb") as f:
        f.write(buf.getvalue())

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.num_points == num_points
    assert len(loaded.extensions) == 1

    ext = loaded.extensions[0]
    assert isinstance(ext, spz.SpzExtensionSafeOrbitCameraAdobe)
    assert abs(ext.safe_orbit_elevation_min - elev_min) < 1e-5
    assert abs(ext.safe_orbit_elevation_max - elev_max) < 1e-5
    assert abs(ext.safe_orbit_radius_min - radius_min) < 1e-5
