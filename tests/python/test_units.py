"""Tests for length units: PackOptions.from_unit, UnpackOptions.to_unit and SPZ_ADOBE_units."""

import math
import struct

import numpy as np
import pytest

import spz

U = spz.LengthUnit
METERS_PER_UNIT = {U.MILLIMETER: 1e-3, U.CENTIMETER: 1e-2, U.METER: 1.0, U.KILOMETER: 1e3}
needs_extensions = pytest.mark.skipif(not spz.has_extension_support(), reason="built without extension support")


def _make_cloud(num_points=20, extent=1.0, scale_range=(-3.0, -1.0), seed=0):
    rng = np.random.default_rng(seed)
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.positions = rng.uniform(-extent, extent, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(*scale_range, size=num_points * 3).astype(np.float32)
    cloud.rotations = np.tile(np.float32([0, 0, 0, 1]), num_points)
    cloud.alphas = np.zeros(num_points, dtype=np.float32)
    cloud.colors = np.zeros(num_points * 3, dtype=np.float32)
    return cloud


def _units(unit):
    ext = spz.SpzExtensionUnitsAdobe()
    ext.length_unit = unit
    return ext


def _pack(from_unit=None):
    options = spz.PackOptions()
    if from_unit is not None:
        options.from_unit = from_unit
    return options


def _unpack(to_unit=None):
    options = spz.UnpackOptions()
    if to_unit is not None:
        options.to_unit = to_unit
    return options


def _of_type(cloud, cls):
    return [e for e in cloud.extensions if isinstance(e, cls)]


def test_options_default_to_meters():
    assert spz.PackOptions().from_unit == U.METER
    assert spz.UnpackOptions().to_unit == U.METER


def test_meter_defaults_are_a_no_op():
    """Without a units record, default options neither change the bytes nor the loaded values."""
    cloud = _make_cloud()
    data = spz.save_spz_to_buffer(cloud, spz.PackOptions())
    assert data == spz.save_spz_to_buffer(cloud, _pack(U.UNSPECIFIED))
    default = spz.load_spz_from_buffer(data)
    raw = spz.load_spz_from_buffer(data, _unpack(U.UNSPECIFIED))
    for attr in ("positions", "scales", "rotations", "alphas", "colors"):
        np.testing.assert_array_equal(getattr(default, attr), getattr(raw, attr))


def test_file_without_record_converts_from_meters():
    """A file with no units record is meters, so to_unit converts from meters."""
    cloud = _make_cloud()
    loaded = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()), _unpack(U.MILLIMETER))
    meters = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    np.testing.assert_allclose(loaded.positions, meters.positions * 1000, rtol=1e-6)
    np.testing.assert_allclose(loaded.scales, meters.scales + math.log(1000), atol=1e-5)


@needs_extensions
@pytest.mark.parametrize("unit", [U.MILLIMETER, U.CENTIMETER, U.METER, U.KILOMETER])
def test_round_trip_in_each_storage_unit(unit):
    """Data stored in any unit loads back in meters by default, and the record survives."""
    cloud = _make_cloud()
    cloud.extensions = [_units(unit)]
    data = spz.save_spz_to_buffer(cloud, spz.PackOptions())
    assert struct.unpack_from("<III", data, 32) == (0xADBE0006, 4, unit.value)

    loaded = spz.load_spz_from_buffer(data)

    grid = METERS_PER_UNIT[unit] / 4096
    np.testing.assert_allclose(loaded.positions, cloud.positions, atol=grid / 2 + 1e-7)
    np.testing.assert_allclose(loaded.scales, cloud.scales, atol=1 / 32 + 1e-5)
    (record,) = _of_type(loaded, spz.SpzExtensionUnitsAdobe)
    assert record.length_unit == unit


@needs_extensions
def test_stored_values_are_in_the_record_unit():
    """Millimeter storage holds positions * 1000 and log scales + ln(1000)."""
    cloud = _make_cloud()
    cloud.extensions = [_units(U.MILLIMETER)]
    stored = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()), _unpack(U.UNSPECIFIED))
    np.testing.assert_allclose(stored.positions, cloud.positions * 1000, atol=0.5 / 4096 + 1e-4)
    np.testing.assert_allclose(stored.scales, cloud.scales + math.log(1000), atol=1 / 32 + 1e-5)


@needs_extensions
@pytest.mark.parametrize("to_unit,factor", [(U.MILLIMETER, 1000.0), (U.CENTIMETER, 100.0), (U.KILOMETER, 1e-3)])
def test_explicit_target_unit(to_unit, factor):
    cloud = _make_cloud()
    cloud.extensions = [_units(U.MILLIMETER)]
    data = spz.save_spz_to_buffer(cloud, spz.PackOptions())
    meters = spz.load_spz_from_buffer(data)
    loaded = spz.load_spz_from_buffer(data, _unpack(to_unit))
    np.testing.assert_allclose(loaded.positions, meters.positions * factor, rtol=1e-6)
    np.testing.assert_allclose(loaded.scales, meters.scales + math.log(factor), atol=1e-5)


@needs_extensions
def test_from_unit_converts_input():
    """A cloud in centimeters with from_unit=CENTIMETER stores and loads the same as its meter twin."""
    meters = _make_cloud()
    centimeters = _make_cloud()
    centimeters.positions = meters.positions * 100
    centimeters.scales = meters.scales + np.float32(math.log(100))
    meters.extensions = [_units(U.MILLIMETER)]
    centimeters.extensions = [_units(U.MILLIMETER)]
    a = spz.load_spz_from_buffer(spz.save_spz_to_buffer(meters, spz.PackOptions()))
    b = spz.load_spz_from_buffer(spz.save_spz_to_buffer(centimeters, _pack(U.CENTIMETER)))
    np.testing.assert_allclose(a.positions, b.positions, atol=1e-6)
    np.testing.assert_array_equal(a.scales, b.scales)


@needs_extensions
def test_millimeters_keep_scales_below_the_meter_clamp():
    """Meters clamp log scales at -10 (45 um); millimeters represent a 20 um splat."""
    cloud = _make_cloud()
    cloud.scales = np.full(cloud.num_points * 3, math.log(20e-6), dtype=np.float32)
    in_meters = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    cloud.extensions = [_units(U.MILLIMETER)]
    in_mm = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))

    np.testing.assert_allclose(in_meters.scales, -10.0)
    np.testing.assert_allclose(in_mm.scales, math.log(20e-6), atol=1 / 32)


@needs_extensions
def test_millimeters_resolve_a_small_object_below_one_micron():
    """For a 0.4 m object, millimeter storage has a 0.24 um grid; meter storage 0.24 mm."""
    cloud = _make_cloud(num_points=2000, extent=0.2, seed=3)
    in_meters = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    cloud.extensions = [_units(U.MILLIMETER)]
    in_mm = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))

    assert np.abs(in_mm.positions - cloud.positions).max() < 1e-6
    assert np.abs(in_meters.positions - cloud.positions).max() > 50e-6


@needs_extensions
def test_safe_orbit_radius_and_environment_origin_convert():
    """Lengths in other extensions follow the units; angles and directions don't."""
    cloud = _make_cloud(4)
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_elevation_min, orbit.safe_orbit_elevation_max, orbit.safe_orbit_radius_min = -0.5, 1.2, 0.3
    pixels = np.zeros(2 * 4 * 3, dtype=np.float32)
    env = spz.SpzExtensionEnvironmentMapAdobe(4, 2, pixels, origin=(0.1, -0.2, 0.3), up=(0.0, 0.6, 0.8),
                                             e1=(1.0, 0.0, 0.0), e2=(0.0, -0.8, 0.6))
    cloud.extensions = [_units(U.MILLIMETER), orbit, env]
    data = spz.save_spz_to_buffer(cloud, spz.PackOptions())

    assert orbit.safe_orbit_radius_min == pytest.approx(0.3)  # saving didn't touch the caller's records
    np.testing.assert_allclose(env.origin, (0.1, -0.2, 0.3))
    for to_unit, factor in [(None, 1.0), (U.UNSPECIFIED, 1000.0), (U.CENTIMETER, 100.0)]:
        loaded = spz.load_spz_from_buffer(data, _unpack(to_unit))
        (loaded_orbit,) = _of_type(loaded, spz.SpzExtensionSafeOrbitCameraAdobe)
        (loaded_env,) = _of_type(loaded, spz.SpzExtensionEnvironmentMapAdobe)
        assert loaded_orbit.safe_orbit_radius_min == pytest.approx(0.3 * factor)
        assert (loaded_orbit.safe_orbit_elevation_min, loaded_orbit.safe_orbit_elevation_max) == pytest.approx((-0.5, 1.2))
        np.testing.assert_allclose(loaded_env.origin, np.array([0.1, -0.2, 0.3]) * factor, rtol=1e-6)
        np.testing.assert_allclose(loaded_env.up, (0.0, 0.6, 0.8), rtol=1e-6)
        np.testing.assert_allclose(loaded_env.e2, (0.0, -0.8, 0.6), rtol=1e-6)


@needs_extensions
def test_spatial_custom_attributes_are_not_scaled():
    """SPATIAL custom attributes are directions, so a unit change leaves them alone."""
    cloud = _make_cloud(4)
    normals = spz.SpzExtensionCustomAttributeNiantic()
    normals.name = "normal"
    normals.dim = 3
    normals.semantic = spz.SpzAttributeSemantic.SPATIAL
    normals.bits = 32
    values = np.float32([[0, 1, 0], [1, 0, 0], [0, 0, 1], [0.6, 0.8, 0]])
    normals.values = values
    cloud.extensions = [_units(U.MILLIMETER), normals]
    options = _unpack(U.MILLIMETER)
    options.custom_attributes = None
    loaded = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()), options)
    np.testing.assert_array_equal(spz.find_custom_attribute(loaded, "normal").values.reshape(4, 3), values)


@needs_extensions
def test_with_coordinate_system_record():
    """Units and coordinate system records combine: RUB meters in, stored RDF millimeters, RUB meters out."""
    cloud = _make_cloud()
    coords = spz.SpzExtensionCoordinateSystemAdobe()
    coords.coordinate_system = spz.RDF
    cloud.extensions = [coords, _units(U.MILLIMETER)]
    pack = _pack()
    pack.from_coord = spz.RUB
    data = spz.save_spz_to_buffer(cloud, pack)

    stored = spz.load_spz_from_buffer(data, _unpack(U.UNSPECIFIED))
    expected = cloud.positions.reshape(-1, 3) * np.float32([1, -1, -1]) * 1000
    np.testing.assert_allclose(stored.positions.reshape(-1, 3), expected, atol=1e-3)

    unpack = _unpack()
    unpack.to_coord = spz.RUB
    np.testing.assert_allclose(spz.load_spz_from_buffer(data, unpack).positions, cloud.positions, atol=1e-6)


@needs_extensions
def test_georeference_applies_to_meter_positions():
    """The georeference transform is defined on RUB meters, so it holds for default loads of any storage unit."""
    from scipy.spatial.transform import Rotation

    cloud = _make_cloud()
    geo = spz.SpzExtensionGeoreferenceNiantic()
    geo.crs = "EPSG:4978"
    geo.origin = [4194304.5, -555555.25, 4713930.125]
    geo.rotation = [0.0, 0.0, 0.70710678118654757, 0.70710678118654757]
    geo.scale = 1.0
    cloud.extensions = [geo, _units(U.MILLIMETER)]
    unpack = _unpack()
    unpack.to_coord = spz.RUB
    loaded = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()), unpack)

    (loaded_geo,) = _of_type(loaded, spz.SpzExtensionGeoreferenceNiantic)
    assert (loaded_geo.origin, loaded_geo.scale) == (geo.origin, geo.scale)  # never scaled
    transform = Rotation.from_quat(loaded_geo.rotation)
    expected = transform.apply(cloud.positions.reshape(-1, 3).astype(np.float64)) + loaded_geo.origin
    actual = transform.apply(loaded.positions.reshape(-1, 3).astype(np.float64)) + loaded_geo.origin
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-6)


@needs_extensions
def test_unspecified_record_means_meters():
    cloud = _make_cloud()
    plain = spz.save_spz_to_buffer(cloud, spz.PackOptions())
    cloud.extensions = [_units(U.UNSPECIFIED)]
    loaded = spz.load_spz_from_buffer(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    np.testing.assert_array_equal(loaded.positions, spz.load_spz_from_buffer(plain).positions)


@needs_extensions
@pytest.mark.parametrize("value", [5, 0xFFFFFFFF])
def test_unknown_unit_record_is_dropped(value):
    """A record with an unknown unit is dropped on load; the data then reads as meters."""
    cloud = _make_cloud()
    cloud.extensions = [_units(U.CENTIMETER)]
    data = bytearray(spz.save_spz_to_buffer(cloud, spz.PackOptions()))
    struct.pack_into("<I", data, 40, value)

    loaded = spz.load_spz_from_buffer(bytes(data))

    assert _of_type(loaded, spz.SpzExtensionUnitsAdobe) == []
    np.testing.assert_allclose(loaded.positions, cloud.positions * 100, atol=1e-3)


def test_ply_is_meters(tmp_path):
    """PLY has no units record: it is written in meters and converted to to_unit on load."""
    cloud = _make_cloud()
    millimeters = _make_cloud()
    millimeters.positions = cloud.positions * 1000
    millimeters.scales = cloud.scales + np.float32(math.log(1000))
    if spz.has_extension_support():
        millimeters.extensions = [_units(U.MILLIMETER)]
    path = str(tmp_path / "units.ply")
    assert spz.save_splat_to_ply(millimeters, _pack(U.MILLIMETER), path)

    in_meters = spz.load_splat_from_ply(path, spz.UnpackOptions())
    np.testing.assert_allclose(in_meters.positions, cloud.positions, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(in_meters.scales, cloud.scales, atol=1e-5)
    if spz.has_extension_support():
        assert in_meters.extensions == []
    in_mm = spz.load_splat_from_ply(path, _unpack(U.MILLIMETER))
    np.testing.assert_allclose(in_mm.positions, millimeters.positions, rtol=1e-6)
