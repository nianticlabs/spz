import numpy as np

from spz.ply_file_ops import gaussian_cloud_to_ply_file, read_ply_to_gaussian_cloud
from spz.spz_file_ops import gaussian_cloud_to_spz_file, read_spz_to_gaussian_cloud
from spz.spz_test_case import SpzTestCase
from spz import GaussianCloud


class TestSpzRead(SpzTestCase):
    def test_read_spz_version2(self):
        gt_file = self.input_test_path() / "spz/shoe_sh3.spz"
        gaussian_cloud_spz = read_spz_to_gaussian_cloud(gt_file)
        num_points = 147856
        self.assertEqual(gaussian_cloud_spz.numPoints, num_points)
        self.assertEqual(len(gaussian_cloud_spz.alphas), num_points)
        self.assertEqual(len(gaussian_cloud_spz.colors), num_points*3)
        self.assertEqual(len(gaussian_cloud_spz.positions), num_points*3)
        self.assertEqual(len(gaussian_cloud_spz.rotations), num_points*4)
        self.assertEqual(len(gaussian_cloud_spz.scales), num_points*3)
        self.assertEqual(len(gaussian_cloud_spz.sh), num_points*15*3)