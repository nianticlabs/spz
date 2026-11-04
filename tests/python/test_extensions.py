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


def _env_map(width=4, height=2, **kwargs):
    pixels = np.linspace(-70.0, 70.0, height * width * 3).astype(np.float16).astype(np.float32)
    pixels[:4] = [0.0, 2.0**-24, 65504.0, -65504.0]  # float16 zero, smallest subnormal, extremes
    args = dict(origin=(0.1, -0.2, 0.3), up=(0.0, 0.6, 0.8), e1=(1.0, 0.0, 0.0), e2=(0.0, -0.8, 0.6),
                metadata='{"format": "deferred_reflection_v1", "env": {"mip_levels": 6}}')
    args.update(kwargs)
    return spz.SpzExtensionEnvironmentMapAdobe(width, height, pixels, **args)


def _env_maps(cloud):
    return [e for e in cloud.extensions if isinstance(e, spz.SpzExtensionEnvironmentMapAdobe)]


def _assert_env_map_equal(loaded, expected):
    assert (loaded.width, loaded.height, loaded.metadata) == (expected.width, expected.height, expected.metadata)
    np.testing.assert_array_equal(loaded.pixels, expected.pixels)
    for name in ("origin", "up", "e1", "e2"):
        np.testing.assert_array_equal(np.float32(getattr(loaded, name)), np.float32(getattr(expected, name)))


def _save_ply(cloud, path, pack_options=None):
    return spz.save_splat_to_ply(cloud, pack_options or spz.PackOptions(), str(path))


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_environment_map_spz_round_trip():
    """Float16 pixels, the frame and the metadata load exactly, next to other extensions."""
    env = _env_map()
    cloud = _make_cloud(20)
    coords = spz.SpzExtensionCoordinateSystemAdobe()
    coords.coordinate_system = spz.RUB
    cloud.extensions = [coords, env, spz.SpzExtensionSafeOrbitCameraAdobe()]

    loaded = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))

    assert loaded.num_points == 20
    assert len(loaded.extensions) == 3
    (loaded_env,) = _env_maps(loaded)
    _assert_env_map_equal(loaded_env, env)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_environment_map_ply_round_trip(tmp_path):
    env = _env_map()
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_radius_min = 0.3
    cloud = _make_cloud(20)
    cloud.extensions = [env, orbit]
    assert _save_ply(cloud, tmp_path / "env.ply")

    loaded = spz.load_splat_from_ply(str(tmp_path / "env.ply"), spz.UnpackOptions())

    assert loaded.num_points == 20
    (loaded_env,) = _env_maps(loaded)
    _assert_env_map_equal(loaded_env, env)
    assert any(isinstance(e, spz.SpzExtensionSafeOrbitCameraAdobe) for e in loaded.extensions)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("fmt,from_coord,to_coord", [
    ("spz", spz.RDF, spz.RUB),
    ("spz", spz.LUF, spz.RDF),
    ("ply", spz.RUB, spz.RUB),
    ("ply", spz.RUB, spz.LUF),
])
def test_environment_map_frame_converts_like_positions(tmp_path, fmt, from_coord, to_coord):
    """origin/up/e1/e2 set to point positions still equal those positions after any conversion."""
    cloud = _make_cloud(4)
    p = cloud.positions.reshape(4, 3)
    env = _env_map(origin=p[0], up=p[1], e1=p[2], e2=p[3])
    cloud.extensions = [env]
    pack_options = spz.PackOptions()
    pack_options.from_coord = from_coord
    unpack_options = spz.UnpackOptions()
    unpack_options.to_coord = to_coord

    if fmt == "spz":
        loaded = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, pack_options), unpack_options)
    else:
        assert _save_ply(cloud, tmp_path / "env.ply", pack_options)
        loaded = spz.load_splat_from_ply(str(tmp_path / "env.ply"), unpack_options)

    (loaded_env,) = _env_maps(loaded)
    frame = np.float32([loaded_env.origin, loaded_env.up, loaded_env.e1, loaded_env.e2])
    np.testing.assert_allclose(frame, loaded.positions.reshape(4, 3), atol=2.0 / 4096)
    if from_coord != to_coord:
        assert not np.allclose(frame, p)
    np.testing.assert_array_equal(np.float32(env.origin), p[0])  # saving didn't touch the caller's record


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_environment_map_ply_with_unknown_elements(tmp_path):
    """Unknown elements before, between and after known ones are skipped by size, in file order."""
    env = _env_map()
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_radius_min = 0.3
    cloud = _make_cloud(20)
    cloud.extensions = [orbit, env]
    assert _save_ply(cloud, tmp_path / "plain.ply")
    data = (tmp_path / "plain.ply").read_bytes()
    header_end = data.index(b"end_header\n")
    header, body = data[:header_end], data[header_end:]
    env_bytes = 56 + env.width * env.height * 12 + len(env.metadata)
    orbit_at = len(body) - env_bytes - 12
    env_at = len(body) - env_bytes

    header = (header
              .replace(b"element safe_orbit", b"element up_axis 1\nproperty uchar up_axis\nelement safe_orbit", 1)
              .replace(b"element environment ", b"element chirality 2\nproperty uint chirality\nelement environment ", 1)
              + b"element trailer 3\nproperty uchar trailer\n")
    body = (body[:orbit_at] + b"\x01" + body[orbit_at:env_at] + struct.pack("<2I", 7, 9) + body[env_at:] + b"xyz")
    (tmp_path / "extra.ply").write_bytes(header + body)

    loaded = spz.load_splat_from_ply(str(tmp_path / "extra.ply"), spz.UnpackOptions())

    assert loaded.num_points == 20
    (loaded_env,) = _env_maps(loaded)
    _assert_env_map_equal(loaded_env, env)
    (loaded_orbit,) = [e for e in loaded.extensions if isinstance(e, spz.SpzExtensionSafeOrbitCameraAdobe)]
    assert loaded_orbit.safe_orbit_radius_min == pytest.approx(0.3)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("make_extensions", [
    lambda: [spz.SpzExtensionEnvironmentMapAdobe(4, 2, np.zeros(4 * 2 * 3 - 1, np.float32))],
    lambda: [spz.SpzExtensionEnvironmentMapAdobe(0, 2, np.zeros(0, np.float32))],
    lambda: [spz.SpzExtensionEnvironmentMapAdobe(1, 1, np.float32([0.0, np.nan, 0.0]))],
    lambda: [spz.SpzExtensionEnvironmentMapAdobe(1, 1, np.float32([0.0, 70000.0, 0.0]))],
    lambda: [spz.SpzExtensionEnvironmentMapAdobe(1, 1, np.zeros(3, np.float32), up=(0.0, np.inf, 0.0))],
    lambda: [_env_map(), _env_map()],
], ids=["pixel-count", "zero-width", "nan-pixel", "pixel-range", "inf-vector", "duplicate"])
@pytest.mark.parametrize("fmt", ["spz", "ply"])
def test_invalid_environment_map_fails_to_save(tmp_path, fmt, make_extensions):
    cloud = _make_cloud(5)
    cloud.extensions = make_extensions()
    if fmt == "spz":
        with pytest.raises(RuntimeError):
            spz.save_spz_to_buffer(cloud, spz.PackOptions())
    else:
        assert not _save_ply(cloud, tmp_path / "bad.ply")


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("offset,value", [
    (40, struct.pack("<I", 5)),            # width: pixels no longer cover width * height
    (44, struct.pack("<I", 0)),            # height 0
    (96, struct.pack("<I", 0xFFFFFFFF)),   # metadataLength past the end of the record
    (60, struct.pack("<f", np.nan)),       # up.x
], ids=["width", "zero-height", "metadata-length", "nan-vector"])
def test_corrupt_environment_map_is_dropped(offset, value):
    """A malformed record is skipped; the core splat data still loads."""
    cloud = _make_cloud(5)
    cloud.extensions = [_env_map()]
    data = bytearray(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    # header (32) + record type and length (8) -> width 40, height 44, origin 48, up 60, metadataLength 96
    data[offset:offset + len(value)] = value

    loaded = spz.load_spz_from_buffer(bytes(data))

    assert loaded.num_points == 5
    assert _env_maps(loaded) == []


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_corrupt_environment_map_ply_is_dropped(tmp_path):
    """A PLY environment map whose pixel count disagrees with its size is dropped; what follows still loads."""
    env = _env_map()
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_radius_min = 0.3
    cloud = _make_cloud(5)
    cloud.extensions = [env, orbit]
    assert _save_ply(cloud, tmp_path / "env.ply")
    data = (tmp_path / "env.ply").read_bytes()
    count = env.width * env.height
    data = data.replace(f"element environment_map {count}\n".encode(), f"element environment_map {count - 1}\n".encode())
    pixels_end = len(data) - 12 - len(env.metadata)
    (tmp_path / "bad.ply").write_bytes(data[:pixels_end - 12] + data[pixels_end:])

    loaded = spz.load_splat_from_ply(str(tmp_path / "bad.ply"), spz.UnpackOptions())

    assert loaded.num_points == 5
    assert _env_maps(loaded) == []
    (loaded_orbit,) = [e for e in loaded.extensions if isinstance(e, spz.SpzExtensionSafeOrbitCameraAdobe)]
    assert loaded_orbit.safe_orbit_radius_min == pytest.approx(0.3)
