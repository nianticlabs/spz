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
        - elevation_min: float (radians)
        - elevation_max: float (radians)
        - radius_min: float
    """
    packed_gaussians = spz.loadSpzPacked(str(filename))  # pylint: disable=no-member
    for ext in packed_gaussians.extensions:
        if isinstance(ext, spz.SpzExtensionSafeOrbitCameraAdobe):
            return {
                "elevation_min": ext.safeOrbitElevationMin,
                "elevation_max": ext.safeOrbitElevationMax,
                "radius_min": ext.safeOrbitRadiusMin,
            }

    return {}


def read_spz_sh_quantization_data(filename: Path) -> dict:
    """Read spherical harmonics quantization data from an SPZ file.

    Returns:
        Dictionary containing SH quantization data with keys:
        - sh1_bits: int
        - sh_rest_bits: int
        - sh_min: float
        - sh_max: float
    """
    packed_gaussians = spz.loadSpzPacked(str(filename))  # pylint: disable=no-member
    for ext in packed_gaussians.extensions:
        if isinstance(ext, spz.SpzExtensionSHQuantizationAdobe):
            return {
                "sh1_bits": ext.sh1Bits,
                "sh_rest_bits": ext.shRestBits,
                "sh_min": ext.shMin,
                "sh_max": ext.shMax,
            }

    return {}


def gaussian_cloud_to_spz_file(
    gaussian_cloud: spz.GaussianCloud,
    filename: Path,
    coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED,
    sh1_bits: int = 5,
    sh_rest_bits: int = 4,
    enable_sh_min_max_scaling: bool = False,
) -> Path:  # pylint: disable=no-member
    filename.parent.mkdir(parents=True, exist_ok=True)
    pack_options = spz.PackOptions()  # pylint: disable=no-member
    setattr(pack_options, "from", coordinate_system)  # 'from' is a Python keyword, so use setattr
    # use the default values for SH bits in version 2.
    pack_options.sh1Bits = sh1_bits
    pack_options.shRestBits = sh_rest_bits
    pack_options.enableSHMinMaxScaling = enable_sh_min_max_scaling
    if pack_options.enableSHMinMaxScaling is False and len(gaussian_cloud.sh) > 0:
        gaussian_cloud.sh = np.array(gaussian_cloud.sh).clip(-1, 1).tolist()  # type: ignore
    spz.saveSpz(gaussian_cloud, pack_options, str(filename))  # pylint: disable=no-member
    return filename


def gaussian_cloud_to_spz_buffer(
    gaussian_cloud: spz.GaussianCloud,
    coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED,
    sh1_bits: int = 5,
    sh_rest_bits: int = 4,
    enable_sh_min_max_scaling: bool = False,
) -> bytes:  # pylint: disable=no-member
    pack_options = spz.PackOptions()  # pylint: disable=no-member
    setattr(pack_options, "from", coordinate_system)  # 'from' is a Python keyword, so use setattr
    # use the default values for SH bits in version 2.
    pack_options.sh1Bits = sh1_bits
    pack_options.shRestBits = sh_rest_bits
    pack_options.enableSHMinMaxScaling = enable_sh_min_max_scaling
    if pack_options.enableSHMinMaxScaling is False and len(gaussian_cloud.sh) > 0:
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

    assert xyz.shape[0] == opacities.shape[0] == scale.shape[0] == rotation.shape[0] == rgb.shape[0], "All arrays must have the same number of points"
    assert scale.shape[1] == 3, "Scale must have shape [N,3]"
    assert rotation.shape[1] == 4, "Rotation must have shape [N,4]"
    assert rgb.shape[1] == 3, "RGB must have shape [N,3]"
    assert sh_degree >= 0, "SH degree must be non-negative"

    norms = np.linalg.norm(rotation, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    rotation = rotation / norms
    return spz.createGaussianCloudFromNumpy(  # type: ignore[attr-defined]
        np.asarray(xyz, dtype=np.float32),
        np.asarray(opacities, dtype=np.float32),
        np.asarray(scale, dtype=np.float32),
        np.asarray(rotation, dtype=np.float32),
        np.asarray(rgb, dtype=np.float32),
        None if f_rest is None else np.asarray(f_rest, dtype=np.float32),
        int(sh_degree),
    )
