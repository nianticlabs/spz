from spz.ply_file_ops import read_ply_to_gaussian_cloud
from spz.spz_test_case import SpzTestCase
from spz import GaussianCloud  # pylint: disable=no-name-in-module


def _assert_gaussian_cloud_correct_size(test_case: SpzTestCase, gaussian_cloud1: GaussianCloud, num_splats: int, sh_degree: int):
    test_case.assertEqual(gaussian_cloud1.numPoints, num_splats)
    test_case.assertEqual(len(gaussian_cloud1.positions), num_splats * 3)
    test_case.assertEqual(len(gaussian_cloud1.scales), num_splats * 3)
    test_case.assertEqual(len(gaussian_cloud1.rotations), num_splats * 4)
    test_case.assertEqual(len(gaussian_cloud1.alphas), num_splats)
    test_case.assertEqual(len(gaussian_cloud1.colors), num_splats * 3)
    test_case.assertEqual(len(gaussian_cloud1.sh), num_splats * ((sh_degree + 1) ** 2 - 1) * 3)
    test_case.assertEqual(gaussian_cloud1.shDegree, sh_degree)


class TestPly(SpzTestCase):
    def test_read_ply_float32_sh3(self):
        gt_file = self.input_test_path() / "ply/shoe_sh3_float32.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        num_splats = 147754
        _assert_gaussian_cloud_correct_size(self, gaussian_cloud, num_splats, 3)

    def test_read_ply_uint16_sh0(self):
        gt_file = self.input_test_path() / "ply/shoe_sh0_uint16.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        num_splats = 162207
        _assert_gaussian_cloud_correct_size(self, gaussian_cloud, num_splats, 0)

    def test_read_ply_uint16_sh3(self):
        gt_file = self.input_test_path() / "ply/tree_sh3_uint16.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        num_splats = 156553
        _assert_gaussian_cloud_correct_size(self, gaussian_cloud, num_splats, 3)

    def test_read_ply_float16_sh3_toy(self):
        gt_file = self.input_test_path() / "ply/Furry_Toy_2_sh3_float16.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        num_splats = 262805
        _assert_gaussian_cloud_correct_size(self, gaussian_cloud, num_splats, 3)

    def test_read_ply_float32_sh3_toy(self):
        gt_file = self.input_test_path() / "ply/Furry_Toy_2_sh3_float32.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        num_splats = 262805
        _assert_gaussian_cloud_correct_size(self, gaussian_cloud, num_splats, 3)

    def test_read_ply_float32_sh4(self):
        gt_file = self.input_test_path() / "ply/plant_sh4_float32.ply"
        gaussian_cloud = read_ply_to_gaussian_cloud(gt_file)
        num_splats = 203800
        _assert_gaussian_cloud_correct_size(self, gaussian_cloud, num_splats, 4)
