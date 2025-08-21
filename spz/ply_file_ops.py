from pathlib import Path

import numpy as np
from plyfile import PlyData

import spz


def _extract_positions(vertex_data, num_points, reinterpret_dtype: str = "float32", cast_dtype: str = "float32"):
    print(f"[SPZ] Extract positions (x, y, z), reinterpret as {reinterpret_dtype}, cast as {cast_dtype}.")
    positions = np.array([(vertex_data[i]["x"], vertex_data[i]["y"], vertex_data[i]["z"]) for i in range(num_points)])
    if reinterpret_dtype != cast_dtype:
        positions = positions.view(reinterpret_dtype).astype(cast_dtype)
    return positions.tolist()


def _extract_scales(vertex_data, num_points, reinterpret_dtype: str = "float32", cast_dtype: str = "float32"):
    print(f"[SPZ] Extract scales (scale_0, scale_1, scale_2), reinterpret as {reinterpret_dtype}, cast as {cast_dtype}.")
    scales = np.array([(vertex_data[i]["scale_0"], vertex_data[i]["scale_1"], vertex_data[i]["scale_2"]) for i in range(num_points)])
    if reinterpret_dtype != cast_dtype:
        scales = scales.view(reinterpret_dtype).astype(cast_dtype)
    return scales.tolist()


def _extract_rotations(vertex_data, num_points, reinterpret_dtype: str = "float32", cast_dtype: str = "float32"):
    print(f"[SPZ] Extract rotations (rot_1, rot_2, rot_3, rot_0), reinterpret as {reinterpret_dtype}, cast as {cast_dtype}.")
    # Extract all rotations as a 2D numpy array (num_points x 4)
    rotations = np.array([(vertex_data[i]["rot_1"], vertex_data[i]["rot_2"], vertex_data[i]["rot_3"], vertex_data[i]["rot_0"]) for i in range(num_points)])
    if reinterpret_dtype != cast_dtype:
        rotations = rotations.view(reinterpret_dtype).astype(cast_dtype)
    # Normalize all quaternions at once using vectorized operations
    norms = np.linalg.norm(rotations, axis=1, keepdims=True)
    normalized_rotations = rotations / norms
    return normalized_rotations.tolist()


def _extract_alphas(vertex_data, num_points, reinterpret_dtype: str = "float32", cast_dtype: str = "float32"):
    print(f"[SPZ] Extract alpha (opacity), reinterpret as {reinterpret_dtype}, cast as {cast_dtype}.")
    alphas = np.array([vertex_data[i]["opacity"] for i in range(num_points)])
    if reinterpret_dtype != cast_dtype:
        alphas = alphas.view(reinterpret_dtype).astype(cast_dtype)
    return alphas.tolist()


def _extract_colors(vertex_data, num_points, reinterpret_dtype: str = "float32", cast_dtype: str = "float32"):
    print(f"[SPZ] Extract colors (f_dc_0, f_dc_1, f_dc_2), reinterpret as {reinterpret_dtype}, cast as {cast_dtype}.")
    colors = np.array([(vertex_data[i]["f_dc_0"], vertex_data[i]["f_dc_1"], vertex_data[i]["f_dc_2"]) for i in range(num_points)])
    if reinterpret_dtype != cast_dtype:
        colors = colors.view(reinterpret_dtype).astype(cast_dtype)
    return colors.tolist()


def _extract_spherical_harmonics(ply_data, vertex_data, num_points, no_sh, reinterpret_dtype: str = "float32", cast_dtype: str = "float32"):
    sh_dim = 0
    sh_degree = 0
    sh_base_indices = []
    spherical_harmonics = []
    if not no_sh:
        print("[SPZ] Checking for spherical harmonics.")
        end_found = False
        i = 0
        while not end_found:
            field_r = f"f_rest_{i}"
            if field_r in ply_data.elements[0]:
                if i % 3 == 0:
                    sh_dim += 1
                    sh_base_indices.append(i // 3)
                if i in (8, 23, 44, 71):
                    sh_degree += 1
            else:
                end_found = True
            i += 1
        if sh_degree > 0:
            print(f"[SPZ] Extract spherical harmonics with degree {sh_degree} and dimension {sh_dim}.")
            for i in range(num_points):
                for j in sh_base_indices:
                    spherical_harmonics.append(vertex_data[i][f"f_rest_{j}"])  # R channel
                    spherical_harmonics.append(vertex_data[i][f"f_rest_{j + sh_dim}"])  # G channel
                    spherical_harmonics.append(vertex_data[i][f"f_rest_{j + 2 * sh_dim}"])  # B channel
    if not no_sh and reinterpret_dtype != cast_dtype:
        spherical_harmonics_array = np.array(spherical_harmonics).view(reinterpret_dtype).astype(cast_dtype)
        spherical_harmonics = spherical_harmonics_array.flatten().tolist()
    return sh_degree, spherical_harmonics


def _extract_optional_metadata(plydata):
    print("[SPZ] Extract metadata.")
    metadata = {}

    def _get_element_by_name(name: str):
        for elem in plydata.elements:
            if elem.name == name:
                return elem
        return None

    # Safe orbit elevation min/max (expecting 2 values)
    elevation_elem = _get_element_by_name("safe_orbit_camera_elevation_min_max_radians")
    if elevation_elem and len(elevation_elem.data) == 2:
        metadata["safe_orbit_camera_elevation_min_max_radians"] = np.array(
            [elevation_elem.data[0][0], elevation_elem.data[1][0]], dtype = np.float32
        )
        print(f"[SPZ] camera elevation min/max: {metadata['safe_orbit_camera_elevation_min_max_radians']}.")

    # Safe orbit radius min (single float)
    radius_elem = _get_element_by_name("safe_orbit_camera_radius_min")
    if radius_elem and len(radius_elem.data) >= 1:
        metadata["safe_orbit_camera_radius_min"] = float(radius_elem.data[0][0])
        print(f"[SPZ] camera radius min: {metadata['safe_orbit_camera_radius_min']}.")

    return metadata


# pylint: disable=too-many-locals
def read_ply_to_gaussian_cloud(ply_filename, no_sh: bool = False) -> spz.GaussianCloud:
    print("[SPZ] Read the PLY file")
    ply_data = PlyData.read(ply_filename)
    vertex_data = ply_data["vertex"]

    print("[SPZ] Extract data into lists")
    num_points = len(vertex_data)

    # Extract attributes
    buffer_to_dtype = {ply_data["vertex"].properties[i].name: ply_data["vertex"].properties[i].dtype()[1:] for i in range(len(ply_data["vertex"].properties))}
    # replace u2 with float16 for compatibility with how flywheel writes the PLY file
    for key, value in buffer_to_dtype.items():
        if value == "u2":
            buffer_to_dtype[key] = "f2"

    positions = _extract_positions(vertex_data, num_points, buffer_to_dtype["x"], "float32")
    scales = _extract_scales(vertex_data, num_points, buffer_to_dtype["scale_0"], "float32")
    rotations = _extract_rotations(vertex_data, num_points, buffer_to_dtype["rot_1"], "float32")
    alphas = _extract_alphas(vertex_data, num_points, buffer_to_dtype["opacity"], "float32")
    colors = _extract_colors(vertex_data, num_points, buffer_to_dtype["f_dc_0"], "float32")
    if "f_rest_0" in buffer_to_dtype:
        sh_degree, spherical_harmonics = _extract_spherical_harmonics(ply_data, vertex_data, num_points, no_sh, buffer_to_dtype["f_rest_0"], "float32")
    else:
        sh_degree = 0
        spherical_harmonics = []

    print("[SPZ] Flatten lists for GaussianCloud")
    positions = [item for sublist in positions for item in sublist]
    scales = [item for sublist in scales for item in sublist]
    rotations = [item for sublist in rotations for item in sublist]
    colors = [item for sublist in colors for item in sublist]

    print("[SPZ] Create a GaussianCloud object")
    gaussian_cloud = spz.GaussianCloud()  # pylint: disable=no-member
    gaussian_cloud.numPoints = num_points
    gaussian_cloud.shDegree = sh_degree
    gaussian_cloud.positions = positions
    gaussian_cloud.scales = scales
    gaussian_cloud.rotations = rotations
    gaussian_cloud.alphas = alphas
    gaussian_cloud.colors = colors
    gaussian_cloud.sh = spherical_harmonics

    metadata = _extract_optional_metadata(ply_data)
    if "safe_orbit_camera_elevation_min_max_radians" in metadata and "safe_orbit_camera_radius_min" in metadata:
        ext = spz.SpzExtensionSafeOrbitCameraAdobe()
        ext.safeOrbitElevationMin = metadata["safe_orbit_camera_elevation_min_max_radians"][0]
        ext.safeOrbitElevationMax = metadata["safe_orbit_camera_elevation_min_max_radians"][1]
        ext.safeOrbitRadiusMin = metadata["safe_orbit_camera_radius_min"]
        gaussian_cloud.extensions = [ext]

    return gaussian_cloud


def gaussian_cloud_to_ply_file(
    gaussian_cloud: spz.GaussianCloud,
    filename: Path,
    coordinate_system: int = spz.CoordinateSystem.UNSPECIFIED,
) -> Path:  # pylint: disable=no-member
    pack_options = spz.PackOptions()  # pylint: disable=no-member
    setattr(pack_options, "from", coordinate_system)  # 'from' is a Python keyword, so use setattr
    spz.saveSplatToPly(gaussian_cloud, pack_options, str(filename))  # pylint: disable=no-member
    return filename
