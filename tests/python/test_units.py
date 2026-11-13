"""Tests for the SPZ_ADOBE_units extension, a label of the length unit the data is written in."""

import struct

import numpy as np
import pytest

import spz

pytestmark = pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")

# A units record as the only extension: [u32 type @32][u32 byteLength @36][u32 unit @40].
_UNIT_OFFSET = 40


def _make_cloud(num_points=20, seed=0):
    rng = np.random.default_rng(seed)
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1
    cloud.positions = rng.uniform(-200.0, 200.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-9.0, 5.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(-1.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * 9).astype(np.float32)
    return cloud


def _units(unit):
    ext = spz.SpzExtensionUnitsAdobe()
    ext.length_unit = unit
    return ext


def _units_records(cloud):
    return [e for e in cloud.extensions if isinstance(e, spz.SpzExtensionUnitsAdobe)]


def _assert_gaussians_equal(a, b):
    assert a.num_points == b.num_points
    for attr in ("positions", "scales", "rotations", "alphas", "colors", "sh"):
        np.testing.assert_array_equal(getattr(a, attr), getattr(b, attr), err_msg=attr)


@pytest.mark.parametrize("name,value", [
    ("UNSPECIFIED", 0), ("MILLIMETER", 1), ("CENTIMETER", 2), ("METER", 3), ("KILOMETER", 4)])
def test_record_round_trips(name, value):
    unit = getattr(spz.LengthUnit, name)
    cloud = _make_cloud()
    cloud.extensions = [_units(unit)]
    data = spz.save_spz_to_buffer(cloud, spz.PackOptions())

    assert struct.unpack_from("<III", data, 32) == (0xADBE0006, 4, value)
    (record,) = _units_records(spz.load_spz_from_buffer(data))
    assert record.length_unit == unit
    assert record.extension_type == spz.SpzExtensionType.SPZ_ADOBE_units


def test_data_loads_exactly_as_stored():
    """The record converts nothing: data loads the same as from a file without it."""
    cloud = _make_cloud()
    plain = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    cloud.extensions = [_units(spz.LengthUnit.MILLIMETER)]
    labeled = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    _assert_gaussians_equal(labeled, plain)


def test_default_record_round_trips_as_unspecified():
    """A record attached without setting a unit saves and loads as UNSPECIFIED."""
    cloud = _make_cloud()
    cloud.extensions = [spz.SpzExtensionUnitsAdobe()]
    data = spz.save_spz_to_buffer(cloud, spz.PackOptions())

    assert struct.unpack_from("<I", data, _UNIT_OFFSET) == (0,)
    (record,) = _units_records(spz.load_spz_from_buffer(data))
    assert record.length_unit == spz.LengthUnit.UNSPECIFIED


def _patched_unit(value):
    """A file holding a units record with the given raw value, followed by a safe orbit record."""
    cloud = _make_cloud()
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_radius_min = 0.3
    cloud.extensions = [_units(spz.LengthUnit.CENTIMETER), orbit]
    data = bytearray(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    baseline = spz.load_spz_from_buffer(bytes(data))
    struct.pack_into("<I", data, _UNIT_OFFSET, value)
    return spz.load_spz_from_buffer(bytes(data)), baseline


def test_unspecified_record_is_kept_on_load():
    loaded, baseline = _patched_unit(0)
    (record,) = _units_records(loaded)
    assert record.length_unit == spz.LengthUnit.UNSPECIFIED
    _assert_gaussians_equal(loaded, baseline)


@pytest.mark.parametrize("value", [5, 0xFFFFFFFF])
def test_unknown_unit_record_is_dropped(value):
    """A record with an unknown unit is dropped on load; the core data and later records still load."""
    loaded, baseline = _patched_unit(value)
    assert _units_records(loaded) == []
    assert [type(e) for e in loaded.extensions] == [spz.SpzExtensionSafeOrbitCameraAdobe]
    _assert_gaussians_equal(loaded, baseline)


@pytest.mark.parametrize("value", [5, 0xFFFFFFFF])
def test_unknown_unit_cannot_be_set(value):
    """Python can't build a record that saveSpz would reject."""
    with pytest.raises(ValueError):
        spz.LengthUnit(value)
    with pytest.raises(TypeError):
        spz.SpzExtensionUnitsAdobe().length_unit = value


def test_ply_drops_the_record(tmp_path):
    cloud = _make_cloud()
    cloud.extensions = [_units(spz.LengthUnit.MILLIMETER)]
    path = str(tmp_path / "units.ply")
    assert spz.save_splat_to_ply(cloud, spz.PackOptions(), path)
    loaded = spz.load_splat_from_ply(path, spz.UnpackOptions())
    assert loaded.extensions == []
    np.testing.assert_array_equal(loaded.positions, cloud.positions)
    np.testing.assert_array_equal(loaded.scales, cloud.scales)
