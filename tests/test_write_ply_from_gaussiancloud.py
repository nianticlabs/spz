import numpy as np

from spz.ply_file_ops import gaussian_cloud_to_ply_file, read_ply_to_gaussian_cloud
from spz.spz_test_case import SpzTestCase
from spz import GaussianCloud  # pylint: disable=no-name-in-module


def _assert_gaussian_cloud_equal(test_case: SpzTestCase, gaussian_cloud1: GaussianCloud, gaussian_cloud2: GaussianCloud):
    test_case.assertEqual(gaussian_cloud1.numPoints, gaussian_cloud2.numPoints)
    test_case.assertEqual(len(gaussian_cloud1.positions), len(gaussian_cloud2.positions))
    test_case.assertEqual(len(gaussian_cloud1.scales), len(gaussian_cloud2.scales))
    test_case.assertEqual(len(gaussian_cloud1.rotations), len(gaussian_cloud2.rotations))
    test_case.assertEqual(len(gaussian_cloud1.alphas), len(gaussian_cloud2.alphas))
    test_case.assertEqual(len(gaussian_cloud1.colors), len(gaussian_cloud2.colors))
    test_case.assertEqual(len(gaussian_cloud1.sh), len(gaussian_cloud2.sh))
    test_case.assertEqual(gaussian_cloud1.shDegree, gaussian_cloud2.shDegree)

    # check the values are close
    test_case.assertTrue(np.allclose(gaussian_cloud1.positions, gaussian_cloud2.positions, atol=1e-3))
    test_case.assertTrue(np.allclose(gaussian_cloud1.scales, gaussian_cloud2.scales, atol=1e-3))
    test_case.assertTrue(np.allclose(gaussian_cloud1.rotations, gaussian_cloud2.rotations, atol=1e-3))
    test_case.assertTrue(np.allclose(gaussian_cloud1.alphas, gaussian_cloud2.alphas, atol=1e-3))
    test_case.assertTrue(np.allclose(gaussian_cloud1.colors, gaussian_cloud2.colors, atol=1e-3))
    test_case.assertTrue(np.allclose(gaussian_cloud1.sh, gaussian_cloud2.sh, atol=1e-3))


class TestPly(NeuralAssetsTestCase):
    def test_write_ply_float32_sh3(self):
        gt_file = self.input_test_path() / "ply/shoe_sh3_float32.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        ply_file = self.output_test_path() / "shoe_sh3_float32.ply"
        gaussian_cloud_to_ply_file(gaussian_cloud, ply_file)
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(ply_file)
        _assert_gaussian_cloud_equal(self, gaussian_cloud, gaussian_cloud_ply)

    def test_write_ply_uint16_sh0(self):
        gt_file = self.input_test_path() / "ply/shoe_sh0_uint16.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        ply_file = self.output_test_path() / "shoe_sh0_uint16.ply"
        gaussian_cloud_to_ply_file(gaussian_cloud, ply_file)
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(ply_file)
        _assert_gaussian_cloud_equal(self, gaussian_cloud, gaussian_cloud_ply)

    def test_write_ply_uint16_sh3(self):
        gt_file = self.input_test_path() / "ply/tree_sh3_uint16.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        ply_file = self.output_test_path() / "shoe_sh0_uint16.ply"
        gaussian_cloud_to_ply_file(gaussian_cloud, ply_file)
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(ply_file)
        _assert_gaussian_cloud_equal(self, gaussian_cloud, gaussian_cloud_ply)

    def test_write_ply_float16_sh3_toy(self):
        gt_file = self.input_test_path() / "ply/Furry_Toy_2_sh3_float16.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        ply_file = self.output_test_path() / "Furry_Toy_2_sh3_float16.ply"
        gaussian_cloud_to_ply_file(gaussian_cloud, ply_file)
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(ply_file)
        _assert_gaussian_cloud_equal(self, gaussian_cloud, gaussian_cloud_ply)

    def test_write_ply_float32_sh3_toy(self):
        gt_file = self.input_test_path() / "ply/Furry_Toy_2_sh3_float32.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        ply_file = self.output_test_path() / "Furry_Toy_2_sh3_float32.ply"
        gaussian_cloud_to_ply_file(gaussian_cloud, ply_file)
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(ply_file)
        _assert_gaussian_cloud_equal(self, gaussian_cloud, gaussian_cloud_ply)

    def test_write_ply_float32_sh4(self):
        gt_file = self.input_test_path() / "ply/plant_sh4_float32.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        ply_file = self.output_test_path() / "plant_sh4_float32.ply"
        gaussian_cloud_to_ply_file(gaussian_cloud, ply_file)
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(ply_file)
        _assert_gaussian_cloud_equal(self, gaussian_cloud, gaussian_cloud_ply)
