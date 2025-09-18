import numpy as np
from test_ply2spz import _assert_gaussian_cloud_equal

from spz.spz_file_ops import gaussian_cloud_from_numpy, gaussian_cloud_to_spz_file, read_spz_to_gaussian_cloud
from spz.spz_test_case import SpzTestCase


class TestPly(SpzTestCase):
    def test_np_to_spz(self):
        n_points = 100000
        xyz = np.random.rand(n_points, 3)
        opacities = np.random.uniform(-2.4, 13.8, n_points) # typical range of opacities
        scale = np.random.uniform(-15.4, -2.5, [n_points, 3]) # typical range of scales
        rotation = np.random.rand(n_points, 4)
        rotation = rotation / np.linalg.norm(rotation, axis=1, keepdims=True)
        rgb = np.random.rand(n_points, 3)
        sh = np.random.rand(n_points, 15*3)

        gaussian_cloud = gaussian_cloud_from_numpy(xyz, opacities, scale, rotation, rgb, sh, 3)
        spz_file = self.output_test_path() / "np_to_spz.spz"
        gaussian_cloud_to_spz_file(gaussian_cloud, spz_file, enable_sh_min_max_scaling=True)
        gaussian_cloud_spz = read_spz_to_gaussian_cloud(spz_file)
        _assert_gaussian_cloud_equal(self, gaussian_cloud, gaussian_cloud_spz)
