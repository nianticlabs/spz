"""Tests for SPZ extension round-trip through NGSP v4 and legacy gzip formats."""

import gzip
import io
import math
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


# ---------------------------------------------------------------------------
# SPZ_NIANTIC_georeference (0x4E530001)
# ---------------------------------------------------------------------------
#
# File layout when the georeference extension is the only extension: fixed 32-byte header,
# then the record at offset 32: [u32 type @32][u32 byteLength @36][payload @40].
# Payload offsets within the file: ext_version @40, flags @41, crs_length @42, origin @44,
# rotation @68, scale @100, epoch @108, crs @116.
_GEOREF_RECORD_OFFSET = 32
_GEOREF_PAYLOAD_OFFSET = 40
_GEOREF_FLAGS_OFFSET = 41
_GEOREF_CRS_LENGTH_OFFSET = 42
_GEOREF_ORIGIN_OFFSET = 44
_GEOREF_EPOCH_OFFSET = 108
_GEOREF_CRS_OFFSET = 116
_GEOREF_MAX_CRS_BYTES = 0xFFFF
# Body-fixed Cartesian CRS for the Moon, which has no geocentric CRS code, given inline as PROJJSON.
_MOON_CARTESIAN_PROJJSON = (
    '{"type":"GeodeticCRS","name":"Moon (2015) - Sphere / Ocentric / Cartesian",'
    '"datum":{"type":"GeodeticReferenceFrame","name":"Moon (2015) - Sphere",'
    '"ellipsoid":{"name":"Moon (2015) - Sphere","radius":1737400},'
    '"prime_meridian":{"name":"Reference Meridian","longitude":0}},'
    '"coordinate_system":{"subtype":"Cartesian","axis":['
    '{"name":"Geocentric X","abbreviation":"X","direction":"geocentricX","unit":"metre"},'
    '{"name":"Geocentric Y","abbreviation":"Y","direction":"geocentricY","unit":"metre"},'
    '{"name":"Geocentric Z","abbreviation":"Z","direction":"geocentricZ","unit":"metre"}]}}'
)


def _make_georef_ext(epoch=float("nan"), scale=0.9736, crs="EPSG:4978", crs_encoding=None):
    ext = spz.SpzExtensionGeoreferenceNiantic()
    ext.crs = crs
    if crs_encoding is not None:
        ext.crs_encoding = crs_encoding
    # Earth-radius magnitudes that are not exactly representable, to check full double precision.
    ext.origin = [6378137.123456789, -6356752.987654321, 1234567.891011121]
    ext.rotation = [0.5, 0.5, 0.5, 0.5]
    ext.scale = scale
    ext.epoch = epoch
    return ext


def _save_to_tmp(cloud, name, pack_options=None):
    filename = os.path.join(tempfile.gettempdir(), name)
    assert spz.save_spz(cloud, pack_options or spz.PackOptions(), filename) is True
    return filename


def _assert_gaussians_equal(a, b):
    assert a.num_points == b.num_points
    for attr in ("positions", "scales", "rotations", "alphas", "colors", "sh"):
        assert np.array_equal(np.asarray(getattr(a, attr)), np.asarray(getattr(b, attr))), (
            f"Gaussian attribute {attr} differs"
        )


def _patch_file(filename, offset, data):
    with open(filename, "r+b") as f:
        f.seek(offset)
        f.write(data)


def _read_flags(filename):
    with open(filename, "rb") as f:
        f.seek(_GEOREF_FLAGS_OFFSET)
        return f.read(1)[0]


def _assert_rejected(filename, baseline):
    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert len(loaded.extensions) == 0
    _assert_gaussians_equal(loaded, baseline)


def _baseline(name):
    return spz.load_spz(_save_to_tmp(_make_cloud(), name), spz.UnpackOptions())


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("epoch", [float("nan"), 2025.5], ids=["no_epoch", "with_epoch"])
def test_georeference_round_trip(epoch):
    """All fields survive a save/load round-trip exactly; reserved flag bits are written as 0."""
    cloud = _make_cloud()
    ext = _make_georef_ext(epoch=epoch)
    cloud.extensions = [ext]

    filename = _save_to_tmp(cloud, "georef_round_trip.spz")
    assert _read_flags(filename) & 0xF8 == 0
    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.num_points == cloud.num_points
    assert len(loaded.extensions) == 1

    e = loaded.extensions[0]
    assert isinstance(e, spz.SpzExtensionGeoreferenceNiantic)
    assert e.crs == ext.crs
    assert e.crs_encoding == spz.CrsEncoding.AUTHORITY_CODE
    assert e.origin == ext.origin
    assert e.rotation == ext.rotation
    assert e.scale == ext.scale
    if math.isnan(epoch):
        assert math.isnan(e.epoch)
    else:
        assert e.epoch == epoch


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_coexists_with_safe_orbit():
    """A file carrying both georeference and safe orbit round-trips both intact."""
    cloud = _make_cloud()
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_elevation_min = -0.5
    orbit.safe_orbit_elevation_max = 1.2
    orbit.safe_orbit_radius_min = 0.3
    georef = _make_georef_ext(epoch=2025.5)
    cloud.extensions = [orbit, georef]

    loaded = spz.load_spz(_save_to_tmp(cloud, "georef_coexist.spz"), spz.UnpackOptions())
    assert len(loaded.extensions) == 2
    loaded_orbit = next(
        e for e in loaded.extensions if isinstance(e, spz.SpzExtensionSafeOrbitCameraAdobe)
    )
    loaded_georef = next(
        e for e in loaded.extensions if isinstance(e, spz.SpzExtensionGeoreferenceNiantic)
    )
    assert abs(loaded_orbit.safe_orbit_elevation_min - orbit.safe_orbit_elevation_min) < 1e-5
    assert abs(loaded_orbit.safe_orbit_elevation_max - orbit.safe_orbit_elevation_max) < 1e-5
    assert abs(loaded_orbit.safe_orbit_radius_min - orbit.safe_orbit_radius_min) < 1e-5
    assert loaded_georef.crs == georef.crs
    assert loaded_georef.origin == georef.origin
    assert loaded_georef.rotation == georef.rotation
    assert loaded_georef.scale == georef.scale
    assert loaded_georef.epoch == georef.epoch


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_coexists_with_non_rub_storage():
    """Storing the data in another frame (RDF) leaves the georeference record unchanged."""
    cloud = _make_cloud()
    coord = spz.SpzExtensionCoordinateSystemAdobe()
    coord.coordinate_system = spz.RDF
    georef = _make_georef_ext(epoch=2025.5)
    cloud.extensions = [coord, georef]
    pack = spz.PackOptions()
    pack.from_coord = spz.CoordinateSystem.RUB

    unpack = spz.UnpackOptions()
    unpack.to_coord = spz.CoordinateSystem.RUB
    loaded = spz.load_spz(_save_to_tmp(cloud, "georef_rdf_storage.spz", pack), unpack)
    e = next(x for x in loaded.extensions if isinstance(x, spz.SpzExtensionGeoreferenceNiantic))
    assert e.crs == georef.crs
    assert e.origin == georef.origin
    assert e.rotation == georef.rotation
    assert e.scale == georef.scale
    assert e.epoch == georef.epoch
    # Positions come back in RUB, the frame the transform is defined from.
    assert np.allclose(np.asarray(loaded.positions), np.asarray(cloud.positions), atol=1e-3)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_skip_path_preserves_core_data():
    """A reader that does not implement the extension still loads identical core data.

    Simulated by rewriting the record's type ID to an unregistered vendor, which
    exercises the same unknown-type skip path a non-implementing reader would take.
    """
    baseline = _baseline("georef_skip_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext()]
    filename = _save_to_tmp(cloud, "georef_skip.spz")
    _patch_file(filename, _GEOREF_RECORD_OFFSET, struct.pack("<I", 0x7A7A0001))
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_truncated_fixed_section_skipped():
    """A record shorter than the fixed section is skipped and core data still loads."""
    baseline = _baseline("georef_trunc_baseline.spz")
    cloud = _make_cloud()
    # A long crs leaves room to carve a filler record out of the same zone.
    cloud.extensions = [
        _make_georef_ext(crs=_MOON_CARTESIAN_PROJJSON, crs_encoding=spz.CrsEncoding.PROJJSON)
    ]
    filename = _save_to_tmp(cloud, "georef_trunc.spz")
    with open(filename, "rb") as f:
        toc_offset = struct.unpack_from("<I", f.read(32), 16)[0]

    # Rewrite the extension zone in place: a georeference record truncated to a 10-byte payload,
    # followed by an unknown-type record padding the remaining zone bytes.
    truncated_len = 10
    filler_offset = _GEOREF_PAYLOAD_OFFSET + truncated_len
    filler_payload_len = toc_offset - filler_offset - 8
    assert filler_payload_len >= 0
    _patch_file(filename, _GEOREF_RECORD_OFFSET + 4, struct.pack("<I", truncated_len))
    _patch_file(filename, filler_offset, struct.pack("<II", 0x7A7A0002, filler_payload_len))
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_crs_length_overrun_skipped():
    """A crs_length exceeding the payload is rejected and core data still loads."""
    baseline = _baseline("georef_crs_overrun_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext()]
    filename = _save_to_tmp(cloud, "georef_crs_overrun.spz")
    _patch_file(filename, _GEOREF_CRS_LENGTH_OFFSET, struct.pack("<H", 0xFFFF))
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_unsupported_version_skipped():
    """A payload with ext_version=2 is skipped and core data still loads."""
    baseline = _baseline("georef_v2_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext()]
    filename = _save_to_tmp(cloud, "georef_v2.spz")
    _patch_file(filename, _GEOREF_PAYLOAD_OFFSET, struct.pack("<B", 2))
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("crs_encoding", [2, 3])
def test_georeference_unknown_crs_encoding_skipped(crs_encoding):
    """An unassigned crs_encoding rejects the record: the target CRS would be unreadable."""
    baseline = _baseline("georef_crs_enc_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext()]
    filename = _save_to_tmp(cloud, "georef_crs_unknown_encoding.spz")
    flags = (_read_flags(filename) & ~0x06) | (crs_encoding << 1)
    _patch_file(filename, _GEOREF_FLAGS_OFFSET, struct.pack("<B", flags))
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_reserved_flag_bits_ignored():
    """Reserved flags bits (3-7) are ignored on read."""
    cloud = _make_cloud()
    ext = _make_georef_ext(epoch=2025.5)
    cloud.extensions = [ext]
    filename = _save_to_tmp(cloud, "georef_reserved_bits.spz")
    _patch_file(filename, _GEOREF_FLAGS_OFFSET, struct.pack("<B", _read_flags(filename) | 0xF8))

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert len(loaded.extensions) == 1
    assert loaded.extensions[0].crs == ext.crs
    assert loaded.extensions[0].epoch == ext.epoch


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize(
    "crs",
    ["", "4978", "EPSG:", ":4978", "EPSG:4978:1", "EPSG 4978", "EPSG:497 8", "EPSG:9988@2026.7"],
    ids=[
        "empty", "no_authority", "no_code", "empty_authority", "two_colons", "space", "inner_space",
        "epoch_suffix",
    ],
)
def test_georeference_malformed_crs_code_skipped(crs):
    """The writer omits a malformed AUTHORITY:CODE, so the record is rejected; core data still loads."""
    baseline = _baseline("georef_crsid_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext(crs=crs)]
    _assert_rejected(_save_to_tmp(cloud, "georef_crsid.spz"), baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize(
    "patched",
    [b"EPSG 4978", b"EPSG:497@", b"EPSG:49:8", b":EPSG4978", b"EPSG4978:", b"EPSG:497\x01", b"EPSG:497\xc3"],
    ids=["space", "at_sign", "two_colons", "empty_authority", "empty_code", "control_byte", "non_ascii"],
)
def test_georeference_reader_rejects_malformed_crs_code(patched):
    """The reader validates the stored code itself, not just what the writer would emit."""
    baseline = _baseline("georef_rd_code_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext(crs="EPSG:4978")]
    filename = _save_to_tmp(cloud, "georef_rd_code.spz")
    assert len(patched) == len("EPSG:4978")
    _patch_file(filename, _GEOREF_CRS_OFFSET, patched)
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize(
    "crs",
    ["EPSG:4978", "EPSG:9988", "IGNF:RGF93"],
    ids=["wgs84_ecef", "itrf2020", "non_numeric"],
)
def test_georeference_crs_code_authorities_round_trip(crs):
    """Geocentric codes from any authority round-trip verbatim — not tied to the EPSG registry."""
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext(crs=crs)]
    loaded = spz.load_spz(_save_to_tmp(cloud, "georef_crsid_authority.spz"), spz.UnpackOptions())
    assert len(loaded.extensions) == 1
    assert loaded.extensions[0].crs == crs


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_projjson_crs_round_trips():
    """A non-Earth target with no registry code round-trips as inline PROJJSON."""
    cloud = _make_cloud()
    cloud.extensions = [
        _make_georef_ext(crs=_MOON_CARTESIAN_PROJJSON, crs_encoding=spz.CrsEncoding.PROJJSON)
    ]
    loaded = spz.load_spz(_save_to_tmp(cloud, "georef_projjson_crs.spz"), spz.UnpackOptions())
    assert len(loaded.extensions) == 1
    assert loaded.extensions[0].crs == _MOON_CARTESIAN_PROJJSON
    assert loaded.extensions[0].crs_encoding == spz.CrsEncoding.PROJJSON


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_crs_at_size_cap_round_trips():
    """A crs of exactly 65,535 bytes (u16 crs_length) round-trips; one byte more is omitted."""
    cloud = _make_cloud()
    at_cap = _MOON_CARTESIAN_PROJJSON + " " * (_GEOREF_MAX_CRS_BYTES - len(_MOON_CARTESIAN_PROJJSON))
    assert len(at_cap) == _GEOREF_MAX_CRS_BYTES

    cloud.extensions = [_make_georef_ext(crs=at_cap, crs_encoding=spz.CrsEncoding.PROJJSON)]
    loaded = spz.load_spz(_save_to_tmp(cloud, "georef_crs_cap.spz"), spz.UnpackOptions())
    assert len(loaded.extensions) == 1
    assert loaded.extensions[0].crs == at_cap

    cloud.extensions = [_make_georef_ext(crs=at_cap + " ", crs_encoding=spz.CrsEncoding.PROJJSON)]
    loaded = spz.load_spz(_save_to_tmp(cloud, "georef_crs_over_cap.spz"), spz.UnpackOptions())
    assert len(loaded.extensions) == 0


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_empty_projjson_crs_skipped():
    """An empty PROJJSON crs is rejected like an empty code; core data still loads."""
    baseline = _baseline("georef_empty_projjson_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext(crs="", crs_encoding=spz.CrsEncoding.PROJJSON)]
    _assert_rejected(_save_to_tmp(cloud, "georef_empty_projjson.spz"), baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize(
    "bad", [b"\xff", b"\xc0\xaf", b"\xed\xa0\x80"], ids=["invalid_byte", "overlong", "surrogate"]
)
def test_georeference_projjson_crs_invalid_utf8_skipped(bad):
    """A PROJJSON crs that is not well-formed UTF-8 rejects the record."""
    baseline = _baseline("georef_projjson_utf8_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [
        _make_georef_ext(crs=_MOON_CARTESIAN_PROJJSON, crs_encoding=spz.CrsEncoding.PROJJSON)
    ]
    filename = _save_to_tmp(cloud, "georef_projjson_bad_utf8.spz")
    _patch_file(filename, _GEOREF_CRS_OFFSET + 1, bad)
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("value", [float("nan"), float("inf")], ids=["nan", "inf"])
def test_georeference_reader_rejects_non_finite_origin(value):
    """A non-finite origin component rejects the record; core data still loads."""
    baseline = _baseline("georef_origin_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext()]
    filename = _save_to_tmp(cloud, "georef_origin_nonfinite.spz")
    _patch_file(filename, _GEOREF_ORIGIN_OFFSET + 8, struct.pack("<d", value))  # origin[1]
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize(
    "rotation",
    [
        [0.0, 0.0, 0.0, 2.0],
        [1e-160, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [float("nan"), 0.0, 0.0, 1.0],
        [float("inf"), 0.0, 0.0, 1.0],
    ],
    ids=["non_unit", "tiny", "zero", "nan", "inf"],
)
def test_georeference_non_unit_rotation_skipped(rotation):
    """A quaternion whose norm is not within 1e-6 of 1 rejects the record; it is never repaired."""
    baseline = _baseline("georef_rot_baseline.spz")
    cloud = _make_cloud()
    ext = _make_georef_ext()
    ext.rotation = rotation
    cloud.extensions = [ext]
    _assert_rejected(_save_to_tmp(cloud, "georef_rot_bad.spz"), baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize(
    "scale", [0.0, -1.0, float("nan"), float("inf")], ids=["zero", "negative", "nan", "inf"]
)
def test_georeference_invalid_scale_skipped(scale):
    """A non-finite or non-positive scale is rejected on load; core data still loads."""
    baseline = _baseline("georef_badscale_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext(scale=scale)]
    _assert_rejected(_save_to_tmp(cloud, "georef_badscale.spz"), baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
@pytest.mark.parametrize("value", [float("nan"), float("inf")], ids=["nan", "inf"])
def test_georeference_reader_rejects_non_finite_present_epoch(value):
    """With has_epoch set, a non-finite epoch rejects the record rather than reading as absent."""
    baseline = _baseline("georef_epoch_baseline.spz")
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext(epoch=2025.5)]
    filename = _save_to_tmp(cloud, "georef_epoch_nonfinite.spz")
    _patch_file(filename, _GEOREF_EPOCH_OFFSET, struct.pack("<d", value))
    _assert_rejected(filename, baseline)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_epoch_absent_when_has_epoch_clear():
    """Without has_epoch, the stored epoch bytes are ignored and epoch reads as NaN."""
    cloud = _make_cloud()
    cloud.extensions = [_make_georef_ext()]
    filename = _save_to_tmp(cloud, "georef_epoch_absent.spz")
    _patch_file(filename, _GEOREF_EPOCH_OFFSET, struct.pack("<d", 2025.5))

    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert len(loaded.extensions) == 1
    assert math.isnan(loaded.extensions[0].epoch)


@pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")
def test_georeference_rejected_record_does_not_affect_following_extension():
    """A rejected georeference record is skipped by length; the next extension still loads."""
    cloud = _make_cloud()
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_radius_min = 0.3
    cloud.extensions = [_make_georef_ext(scale=-1.0), orbit]

    loaded = spz.load_spz(_save_to_tmp(cloud, "georef_bad_then_orbit.spz"), spz.UnpackOptions())
    assert len(loaded.extensions) == 1
    assert isinstance(loaded.extensions[0], spz.SpzExtensionSafeOrbitCameraAdobe)
    assert abs(loaded.extensions[0].safe_orbit_radius_min - 0.3) < 1e-5
