"""Tests for coordinate system conversions and transformations."""

import os
import tempfile

import numpy as np
import pytest

import spz
from test_utils import make_single_point_cloud


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


def test_cross_family_position_rotation():
    """Cross-family convert_coordinates (RUB→RBD) applies 90-degree x-axis rotation."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    # Identity quaternion (xyzw): no rotation.
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)

    # RUB→RBD: flip y,z (→ (1,-2,-3)) then rotate (x,-z,y) (→ (1,3,-2)).
    cloud.convert_coordinates(spz.RUB, spz.RBD)
    np.testing.assert_allclose(
        cloud.positions,
        np.array([1.0, 3.0, -2.0], dtype=np.float32),
        atol=1e-6,
    )

    # Identity quaternion with zero xyz is unaffected by the flip; left-multiplying
    # by the 90-deg-about-x quaternion [s,0,0,s] (xyzw) gives [s,0,0,s].
    s = np.float32(np.sqrt(2.0) / 2.0)
    np.testing.assert_allclose(
        cloud.rotations,
        np.array([s, 0.0, 0.0, s], dtype=np.float32),
        atol=1e-6,
    )


@pytest.mark.parametrize("from_coord,to_coord", [
    # standard→rotated: non-trivial inner flip (RFU↔RBD have flipP=(1,-1,-1))
    (spz.RUB, spz.RBD),
    # standard→rotated: identity inner flip (RUB and RFU are a direct mapping pair)
    (spz.RUB, spz.RFU),
    # left-handed variant of the non-trivial-flip case
    (spz.LUB, spz.LBD),
    # cross-family AND cross-handedness: inner has flipP=(-1,-1,-1), exercises the x-flip path
    (spz.RUB, spz.LBD),
    # rotated→standard (backward direction of each structural case above)
    (spz.RBD, spz.RUB),
    (spz.RFU, spz.RUB),
    (spz.LBD, spz.RUB),
])
@pytest.mark.parametrize("sh_degree", range(5))
def test_cross_family_round_trip(from_coord, to_coord, sh_degree):
    """All fields are restored after a cross-family A→B→A round-trip."""
    rng = np.random.default_rng(seed=42)

    num_sh_coeffs = sum(2 * l + 1 for l in range(1, sh_degree + 1)) * 3

    cloud = spz.GaussianCloud()
    cloud.sh_degree = sh_degree
    cloud.positions = rng.random(3).astype(np.float32)
    q = rng.random(4).astype(np.float32)
    cloud.rotations = q / np.linalg.norm(q)
    cloud.sh = rng.random(num_sh_coeffs).astype(np.float32) if num_sh_coeffs else np.array([], dtype=np.float32)

    original_positions = cloud.positions.copy()
    original_rotations = cloud.rotations.copy()
    original_sh = cloud.sh.copy()

    cloud.convert_coordinates(from_coord, to_coord)
    cloud.convert_coordinates(to_coord, from_coord)

    np.testing.assert_allclose(cloud.positions, original_positions, atol=1e-5)
    np.testing.assert_allclose(cloud.rotations, original_rotations, atol=1e-5)
    if num_sh_coeffs:
        np.testing.assert_allclose(cloud.sh, original_sh, atol=1e-5)


def test_io_rotated_from_coord():
    """Pack with a rotated from_coord; loaded data is in the internal RUB frame."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.antialiased = False
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.sh = np.array([], dtype=np.float32)

    filename = os.path.join(tempfile.gettempdir(), "rotated_from_coord_test.spz")

    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RBD
    assert spz.save_spz(cloud, pack_opts, filename) is True

    # Load with UNSPECIFIED → no output conversion; data is in RUB.
    # Pack converted RBD→RUB: no flip (identity inner), rotate (x,-z,y): (1,2,3)→(1,-3,2).
    loaded = spz.load_spz(filename, spz.UnpackOptions())
    assert loaded.num_points == cloud.num_points
    np.testing.assert_allclose(
        loaded.positions,
        np.array([1.0, -3.0, 2.0], dtype=np.float32),
        atol=1 / 2048.0,
    )


def test_io_rotated_to_coord():
    """Pack with RUB, load with a rotated to_coord; verify position conversion."""
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.antialiased = False
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.sh = np.array([], dtype=np.float32)

    filename = os.path.join(tempfile.gettempdir(), "rotated_to_coord_test.spz")

    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RUB
    assert spz.save_spz(cloud, pack_opts, filename) is True

    # Load with to_coord=RBD → unpack converts RUB→RBD on stored (1,2,3):
    # flip {+1,-1,-1}: (1,-2,-3), rotate (x,-z,y): (1,3,-2).
    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RBD
    loaded = spz.load_spz(filename, unpack_opts)
    assert loaded.num_points == cloud.num_points
    np.testing.assert_allclose(
        loaded.positions,
        np.array([1.0, 3.0, -2.0], dtype=np.float32),
        atol=1 / 2048.0,
    )


def test_cross_family_quaternion_nontrivial():
    """Quaternion rotation is correct for cross-family conversions with a non-trivial quaternion."""
    from scipy.spatial.transform import Rotation as R

    rng = np.random.default_rng(seed=123)
    r_rx_plus = R.from_euler('x', 90, degrees=True)

    # Case 1: RUF→RBU — identity inner flip (flipQ={1,1,1}), pure R_x(+π/2).
    q = rng.random(4).astype(np.float32)
    q /= np.linalg.norm(q)
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.positions = np.array([1.0, 0.0, 0.0], dtype=np.float32)
    cloud.rotations = q.copy()
    cloud.sh = np.array([], dtype=np.float32)
    cloud.convert_coordinates(spz.RUF, spz.RBU)

    expected = (r_rx_plus * R.from_quat(q.astype(np.float64))).as_quat()
    result = cloud.rotations.astype(np.float64)
    if np.dot(result, expected) < 0:
        expected = -expected
    np.testing.assert_allclose(result, expected, atol=1e-5)

    # Case 2: RUB→RBD — R_x(+π/2) first, then inner flip flipQ={1,-1,-1} (from within-rotated RFU↔RBD).
    # Order: rotate-then-flip (forward direction).
    q2 = rng.random(4).astype(np.float32)
    q2 /= np.linalg.norm(q2)
    cloud2 = spz.GaussianCloud()
    cloud2.sh_degree = 0
    cloud2.positions = np.array([1.0, 0.0, 0.0], dtype=np.float32)
    cloud2.rotations = q2.copy()
    cloud2.sh = np.array([], dtype=np.float32)
    cloud2.convert_coordinates(spz.RUB, spz.RBD)

    # Rotate first, then apply flipQ=(1,-1,-1): negate y and z quaternion components.
    q2_rotated = (r_rx_plus * R.from_quat(q2.astype(np.float64))).as_quat()
    expected2 = np.array([q2_rotated[0], -q2_rotated[1], -q2_rotated[2], q2_rotated[3]], dtype=np.float64)
    result2 = cloud2.rotations.astype(np.float64)
    if np.dot(result2, expected2) < 0:
        expected2 = -expected2
    np.testing.assert_allclose(result2, expected2, atol=1e-5)


def test_coordinate_system_extension_overrides_pack_to():
    """Extension on the cloud overrides the 'to' coordinate of packGaussians (default RUB → ext coord)."""
    if not spz.has_extension_support():
        return
    # Cloud positions are in RUB. Attach extension requesting storage in RDF.
    cloud = make_single_point_cloud()
    coord_ext = spz.SpzExtensionCoordinateSystemAdobe()
    coord_ext.coordinate_system = spz.RDF
    cloud.extensions = [coord_ext]

    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RUB  # input is RUB; extension says store as RDF

    filename = os.path.join(tempfile.gettempdir(), "coord_ext_override_test.spz")
    assert spz.save_spz(cloud, pack_opts, filename) is True
    # Load back with no output conversion: extension drives from=RDF, to=UNSPECIFIED → data in RDF.
    loaded = spz.load_spz(filename, spz.UnpackOptions())

    # RUB→RDF flips Y and Z: (1,2,3) → (1,-2,-3).
    np.testing.assert_allclose(
        loaded.positions,
        np.array([1.0, -2.0, -3.0], dtype=np.float32),
        atol=1 / 2048.0,
    )
    # Extension must be preserved in the loaded cloud.
    coord_exts = [
        e for e in loaded.extensions
        if e.extension_type == spz.SpzExtensionType.SPZ_ADOBE_coordinate_system
    ]
    assert len(coord_exts) == 1
    assert coord_exts[0].coordinate_system == spz.RDF


def test_coordinate_system_extension_absent_uses_rub():
    """Without the extension, packGaussians converts to RUB as usual."""
    if not spz.has_extension_support():
        return
    cloud = make_single_point_cloud()
    # No coordinate system extension on the cloud.

    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RDF  # input is RDF; no extension → converts to RUB

    filename = os.path.join(tempfile.gettempdir(), "coord_ext_absent_test.spz")
    assert spz.save_spz(cloud, pack_opts, filename) is True
    # Load back with no output conversion: no extension → from=RUB, to=UNSPECIFIED → data in RUB.
    loaded = spz.load_spz(filename, spz.UnpackOptions())

    # RDF→RUB flips Y and Z: (1,2,3) → (1,-2,-3).
    np.testing.assert_allclose(
        loaded.positions,
        np.array([1.0, -2.0, -3.0], dtype=np.float32),
        atol=1 / 2048.0,
    )


def test_coordinate_system_extension_used_on_unpack():
    """Extension's coordinate is used as 'from' during unpack instead of RUB."""
    if not spz.has_extension_support():
        return
    # Store cloud in RDF (input RUB, extension says store as RDF).
    cloud = make_single_point_cloud()
    coord_ext = spz.SpzExtensionCoordinateSystemAdobe()
    coord_ext.coordinate_system = spz.RDF
    cloud.extensions = [coord_ext]

    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RUB

    filename = os.path.join(tempfile.gettempdir(), "coord_ext_unpack_test.spz")
    assert spz.save_spz(cloud, pack_opts, filename) is True

    # Load with to_coord=RUB: extension drives from=RDF, RDF→RUB flips Y and Z back.
    # Stored value is (1,-2,-3); converting RDF→RUB gives (1,2,3) — original positions.
    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RUB
    loaded = spz.load_spz(filename, unpack_opts)
    np.testing.assert_allclose(
        loaded.positions,
        np.array([1.0, 2.0, 3.0], dtype=np.float32),
        atol=1 / 2048.0,
    )


def test_coordinate_system_extension_unspecified_uses_rub():
    """Extension with UNSPECIFIED coordinateSystem is ignored; packGaussians falls back to RUB."""
    if not spz.has_extension_support():
        return
    cloud = make_single_point_cloud()
    coord_ext = spz.SpzExtensionCoordinateSystemAdobe()
    coord_ext.coordinate_system = spz.UNSPECIFIED  # present but UNSPECIFIED → treated as no extension
    cloud.extensions = [coord_ext]

    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RDF  # input is RDF; extension is UNSPECIFIED → converts to RUB

    filename = os.path.join(tempfile.gettempdir(), "coord_ext_unspecified_test.spz")
    assert spz.save_spz(cloud, pack_opts, filename) is True
    loaded = spz.load_spz(filename, spz.UnpackOptions())

    # Same result as no extension: RDF→RUB flips Y and Z: (1,2,3) → (1,-2,-3).
    np.testing.assert_allclose(
        loaded.positions,
        np.array([1.0, -2.0, -3.0], dtype=np.float32),
        atol=1 / 2048.0,
    )


def test_coordinate_system_extension_rotated_family():
    """Extension with a rotated-family coordinate exercises the cross-family pack/unpack path."""
    if not spz.has_extension_support():
        return
    # Input is in RUB; extension requests storage in RBD (rotated family → cross-family conversion).
    cloud = make_single_point_cloud()
    coord_ext = spz.SpzExtensionCoordinateSystemAdobe()
    coord_ext.coordinate_system = spz.RBD
    cloud.extensions = [coord_ext]

    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RUB

    filename = os.path.join(tempfile.gettempdir(), "coord_ext_rotated_test.spz")
    assert spz.save_spz(cloud, pack_opts, filename) is True

    # Load with to_coord=RUB: extension drives from=RBD, cross-family RBD→RUB restores originals.
    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RUB
    loaded = spz.load_spz(filename, unpack_opts)
    np.testing.assert_allclose(
        loaded.positions,
        np.array([1.0, 2.0, 3.0], dtype=np.float32),
        atol=1 / 2048.0,
    )


def test_cross_family_multipoint():
    """Cross-family convert_coordinates applies the same transform to every point."""
    N = 4
    rng = np.random.default_rng(seed=99)
    cloud = spz.GaussianCloud()
    cloud.sh_degree = 0
    cloud.positions = rng.random(N * 3).astype(np.float32)
    q = rng.random((N, 4)).astype(np.float32)
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    cloud.rotations = q.flatten()
    cloud.sh = np.array([], dtype=np.float32)

    original_positions = cloud.positions.copy()
    original_rotations = cloud.rotations.copy()

    cloud.convert_coordinates(spz.RUB, spz.RBD)
    cloud.convert_coordinates(spz.RBD, spz.RUB)

    np.testing.assert_allclose(cloud.positions, original_positions, atol=1e-5)
    np.testing.assert_allclose(cloud.rotations, original_rotations, atol=1e-5)


def test_ply_cross_family_round_trip():
    """saveSplatToPly/loadSplatFromPly with a cross-family coordinate restores all fields."""
    rng = np.random.default_rng(seed=7)
    N = 3
    sh_degree = 1
    num_sh_per_point = sum(2 * l + 1 for l in range(1, sh_degree + 1)) * 3  # 9

    cloud = spz.GaussianCloud()
    cloud.sh_degree = sh_degree
    cloud.antialiased = False
    cloud.positions = rng.random(N * 3).astype(np.float32)
    cloud.scales = rng.random(N * 3).astype(np.float32)
    q = rng.random((N, 4)).astype(np.float32)
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    cloud.rotations = q.flatten()
    cloud.alphas = rng.random(N).astype(np.float32)
    cloud.colors = rng.random(N * 3).astype(np.float32)
    cloud.sh = rng.random(N * num_sh_per_point).astype(np.float32)

    filename = os.path.join(tempfile.gettempdir(), "ply_cross_family_test.ply")
    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RBD
    assert spz.save_splat_to_ply(cloud, pack_opts, filename) is True

    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RBD
    loaded = spz.load_splat_from_ply(filename, unpack_opts)

    assert loaded.num_points == N
    np.testing.assert_allclose(loaded.positions, cloud.positions, atol=1e-5)
    np.testing.assert_allclose(loaded.sh, cloud.sh, atol=1e-5)
    # PLY stores quaternions as float32 with no quantization; check per-point
    # after resolving the double-cover sign ambiguity.
    for i in range(N):
        r = loaded.rotations[i * 4:(i + 1) * 4]
        e = cloud.rotations[i * 4:(i + 1) * 4]
        if np.dot(r, e) < 0:
            r = -r
        np.testing.assert_allclose(r, e, atol=1e-5)


def test_spz_cross_family_rotation_and_sh():
    """SPZ pack/unpack with a cross-family coordinate correctly transforms rotations and SH.

    Positions are already verified by test_io_rotated_from_coord / test_io_rotated_to_coord.
    Multi-point stride is covered by test_cross_family_multipoint.
    """
    rng = np.random.default_rng(seed=11)
    q = rng.random(4).astype(np.float32)
    q /= np.linalg.norm(q)
    sh = rng.random(9).astype(np.float32)  # sh_degree=1: 3 coeffs × 3 channels

    cloud = spz.GaussianCloud()
    cloud.sh_degree = 1
    cloud.antialiased = False
    cloud.positions = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    cloud.scales = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.rotations = q.copy()
    cloud.alphas = np.array([0.5], dtype=np.float32)
    cloud.colors = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    cloud.sh = sh.copy()

    filename = os.path.join(tempfile.gettempdir(), "spz_cross_family_rot_sh.spz")
    pack_opts = spz.PackOptions()
    pack_opts.from_coord = spz.RBD
    assert spz.save_spz(cloud, pack_opts, filename) is True

    unpack_opts = spz.UnpackOptions()
    unpack_opts.to_coord = spz.RBD
    loaded = spz.load_spz(filename, unpack_opts)

    # SH degree-1 uses 5-bit quantization → max error ~1/16.
    np.testing.assert_allclose(loaded.sh, sh, atol=0.1)
    # Quaternion SmallestThree: 8-bit per component → max error ~1/128.
    r = loaded.rotations
    if np.dot(r, q) < 0:
        r = -r
    np.testing.assert_allclose(r, q, atol=1e-2)
