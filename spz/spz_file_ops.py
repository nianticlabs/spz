from pathlib import Path

import numpy as np
import numpy.typing as npt

import spz


def read_spz_to_gaussian_cloud(filename: Path, coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED) -> spz.GaussianCloud:  # pylint: disable=no-member
    unpack_options = spz.UnpackOptions()  # pylint: disable=no-member
    unpack_options.to = coordinate_system
    gaussian_cloud = spz.loadSpz(str(filename), unpack_options)  # pylint: disable=no-member
    return gaussian_cloud


def read_spz_safe_orbit_data(filename: Path) -> dict:
    """Read safe orbit camera data from an SPZ file.

    Returns:
        Dictionary containing safe orbit data with keys:
        - has_safe_orbit: bool
        - elevation_min: float (radians)
        - elevation_max: float (radians)
        - radius_min: float
    """
    packed_gaussians = spz.loadSpzPacked(str(filename))  # pylint: disable=no-member
    return {
        "has_safe_orbit": packed_gaussians.hasSafeOrbit,
        "elevation_min": packed_gaussians.safeOrbitElevationMin,
        "elevation_max": packed_gaussians.safeOrbitElevationMax,
        "radius_min": packed_gaussians.safeOrbitRadiusMin,
    }


def gaussian_cloud_to_spz_file(
    gaussian_cloud: spz.GaussianCloud,
    filename: Path,
    version: int = spz.LATEST_SPZ_HEADER_VERSION,  # pylint: disable=no-member
    coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED,
    sh1_bits: int = 5,
    sh_rest_bits: int = 5,
    disable_sh_min_max_scaling: bool = False,
    has_safe_orbit: bool = False,
    safe_orbit_elevation_min: float = 0.0,
    safe_orbit_elevation_max: float = 0.0,
    safe_orbit_radius_min: float = 0.0,
) -> Path:  # pylint: disable=no-member
    pack_options = spz.PackOptions()  # pylint: disable=no-member
    pack_options.version = version
    setattr(pack_options, "from", coordinate_system)  # 'from' is a Python keyword, so use setattr
    # use the default values for SH bits in version 2.
    pack_options.sh1Bits = 5 if version < 3 else sh1_bits
    pack_options.shRestBits = 4 if version < 3 else sh_rest_bits
    pack_options.hasSafeOrbit = False if version < 4 else has_safe_orbit
    pack_options.disableSHMinMaxScaling = True if version < 3 else disable_sh_min_max_scaling
    pack_options.safeOrbitElevationMin = safe_orbit_elevation_min
    pack_options.safeOrbitElevationMax = safe_orbit_elevation_max
    pack_options.safeOrbitRadiusMin = safe_orbit_radius_min
    if pack_options.disableSHMinMaxScaling is True and len(gaussian_cloud.sh) > 0:
        gaussian_cloud.sh = np.array(gaussian_cloud.sh).clip(-1, 1).tolist()  # type: ignore
    spz.saveSpz(gaussian_cloud, pack_options, str(filename))  # pylint: disable=no-member
    return filename


def gaussian_cloud_to_spz_buffer(
    gaussian_cloud: spz.GaussianCloud,
    version: int = spz.LATEST_SPZ_HEADER_VERSION,
    coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED,
    sh1_bits: int = 5,
    sh_rest_bits: int = 5,
    disable_sh_min_max_scaling: bool = False,
    has_safe_orbit: bool = False,
    safe_orbit_elevation_min: float = 0.0,
    safe_orbit_elevation_max: float = 0.0,
    safe_orbit_radius_min: float = 0.0,
) -> bytes:  # pylint: disable=no-member
    pack_options = spz.PackOptions()  # pylint: disable=no-member
    pack_options.version = version
    setattr(pack_options, "from", coordinate_system)  # 'from' is a Python keyword, so use setattr
    # use the default values for SH bits in version 2.
    pack_options.sh1Bits = 5 if version < 3 else sh1_bits
    pack_options.shRestBits = 4 if version < 3 else sh_rest_bits
    pack_options.hasSafeOrbit = False if version < 4 else has_safe_orbit
    pack_options.disableSHMinMaxScaling = True if version < 3 else disable_sh_min_max_scaling
    pack_options.safeOrbitElevationMin = safe_orbit_elevation_min
    pack_options.safeOrbitElevationMax = safe_orbit_elevation_max
    pack_options.safeOrbitRadiusMin = safe_orbit_radius_min
    if (version < 3 or disable_sh_min_max_scaling) and len(gaussian_cloud.sh) > 0:
        gaussian_cloud.sh = np.array(gaussian_cloud.sh).clip(-1, 1).tolist()  # type: ignore
    return spz.saveSpzToBytes(gaussian_cloud, pack_options)  # pylint: disable=no-member


def gaussian_cloud_from_numpy(
    xyz: npt.NDArray[np.float32],
    opacities: npt.NDArray[np.float32],
    scale: npt.NDArray[np.float32],
    rotation: npt.NDArray[np.float32],
    rgb: npt.NDArray[np.float32],
    f_rest: npt.NDArray[np.float32] | None = None,
    sh_degree: int = 0,
) -> spz.GaussianCloud:  # pylint: disable=no-member

    gaussian_cloud = spz.GaussianCloud()  # pylint: disable=no-member
    gaussian_cloud.numPoints = xyz.shape[0]
    gaussian_cloud.positions = xyz.flatten().astype(np.float32).tolist()
    gaussian_cloud.colors = rgb.flatten().astype(np.float32).tolist()
    gaussian_cloud.alphas = opacities.flatten().astype(np.float32).tolist()
    gaussian_cloud.scales = scale.flatten().astype(np.float32).tolist()
    # normalize rotations
    print("[SPZ] Normalizing rotations")
    rotation = rotation / np.linalg.norm(rotation, axis=1, keepdims=True)
    gaussian_cloud.rotations = rotation.flatten().astype(np.float32).tolist()
    gaussian_cloud.sh = f_rest.flatten().astype(np.float32).tolist() if f_rest is not None else []
    gaussian_cloud.shDegree = sh_degree

    return gaussian_cloud
