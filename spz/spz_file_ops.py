from pathlib import Path

import numpy as np
import numpy.typing as npt

import spz # pylint: disable=consider-using-from-import


def read_spz_to_gaussian_cloud(filename: Path, coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED) -> spz.GaussianCloud:  # pylint: disable=no-member
    unpack_options = spz.UnpackOptions()  # pylint: disable=no-member
    unpack_options.to = coordinate_system
    gaussian_cloud = spz.loadSpz(str(filename), unpack_options)  # pylint: disable=no-member
    return gaussian_cloud


def gaussian_cloud_to_spz_file(
    gaussian_cloud: spz.GaussianCloud,
    filename: Path,
    coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED,
    sh1_bits: int = 5,
    sh_rest_bits: int = 5,
    disable_sh_min_max_scaling: bool = False,
) -> Path:  # pylint: disable=no-member
    pack_options = spz.PackOptions()  # pylint: disable=no-member
    setattr(pack_options, "from", coordinate_system)  # 'from' is a Python keyword, so use setattr
    setattr(pack_options, "sh1Bits", sh1_bits)
    setattr(pack_options, "shRestBits", sh_rest_bits)
    if disable_sh_min_max_scaling and len(gaussian_cloud.sh) > 0:
        # artificially set the first two sh coefficients to -1 and 1 so the minmax scaler doesn't do anything
        gaussian_cloud.sh = np.array(gaussian_cloud.sh).clip(-1, 1).tolist()
        if np.min(gaussian_cloud.sh) > -1.0:
            gaussian_cloud.sh[0] = -1.0
        if np.max(gaussian_cloud.sh) < 1.0:
            gaussian_cloud.sh[1] = 1.0
    spz.saveSpz(gaussian_cloud, pack_options, str(filename))  # pylint: disable=no-member
    return filename


def gaussian_cloud_to_spz_buffer(
    gaussian_cloud: spz.GaussianCloud,
    coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED,
    sh1_bits: int = 5,
    sh_rest_bits: int = 5,
    disable_sh_min_max_scaling: bool = False,
) -> bytes:  # pylint: disable=no-member
    pack_options = spz.PackOptions()  # pylint: disable=no-member
    setattr(pack_options, "from", coordinate_system)  # 'from' is a Python keyword, so use setattr
    setattr(pack_options, "sh1Bits", sh1_bits)
    setattr(pack_options, "shRestBits", sh_rest_bits)
    if disable_sh_min_max_scaling and len(gaussian_cloud.sh) > 0:
        gaussian_cloud.sh = np.array(gaussian_cloud.sh).clip(-1, 1).tolist()
        # artificially set the first two sh coefficients to -1 and 1 so the minmax scaler doesn't do anything
        if np.min(gaussian_cloud.sh) > -1.0:
            gaussian_cloud.sh[0] = -1.0
        if np.max(gaussian_cloud.sh) < 1.0:
            gaussian_cloud.sh[1] = 1.0
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
