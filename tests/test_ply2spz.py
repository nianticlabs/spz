import numpy as np

from spz import GaussianCloud  # pylint: disable=no-name-in-module
from spz.ply_file_ops import gaussian_cloud_to_ply_file, read_ply_to_gaussian_cloud
from spz.spz_file_ops import gaussian_cloud_to_spz_file, read_spz_to_gaussian_cloud
from spz.spz_test_case import SpzTestCase

COLOR_SCALE = 0.15

MAX_POSITION_ERROR = 0.001
MAX_SCALE_ERROR = 0.015
MAX_ROTATION_ERROR = 0.0075
MAX_ALPHA_ERROR = 0.0025
MAX_COLOR_ERROR = 1.0  # This is in [0,255] space
MAX_SH_ERROR = 0.067


def sigmoid(input_array: np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(-input_array))


# pylint: disable=too-many-locals
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
    positions1 = np.array(gaussian_cloud1.positions)
    positions2 = np.array(gaussian_cloud2.positions)
    pos_diff = np.abs(positions1 - positions2)
    max_pos_diff_idx = np.argmax(pos_diff)
    test_case.assertTrue(
        np.allclose(gaussian_cloud1.positions, gaussian_cloud2.positions, atol=MAX_POSITION_ERROR),
        f"Max position diff {pos_diff[max_pos_diff_idx]} at index {max_pos_diff_idx}: {positions1[max_pos_diff_idx]} vs {positions2[max_pos_diff_idx]}",
    )

    scales1 = np.exp(np.array(gaussian_cloud1.scales))
    scales2 = np.exp(np.array(gaussian_cloud2.scales))
    scale_diff = np.abs(scales1 - scales2)
    max_scale_diff_idx = np.argmax(scale_diff)
    test_case.assertTrue(
        np.allclose(scales1, scales2, atol=MAX_SCALE_ERROR),
        f"Max scale diff {scale_diff[max_scale_diff_idx]} at index {max_scale_diff_idx}: {scales1[max_scale_diff_idx]} vs {scales2[max_scale_diff_idx]}",
    )

    q1s = np.asarray(gaussian_cloud1.rotations).reshape(-1, 4)
    q2s = np.asarray(gaussian_cloud2.rotations).reshape(-1, 4)
    q1s = q1s / np.linalg.norm(q1s, axis=1, keepdims=True)
    q2s = q2s / np.linalg.norm(q2s, axis=1, keepdims=True)
    rotations_dot = np.abs(np.sum(q1s * q2s, axis=1))
    rot_diff = np.abs(rotations_dot - np.ones(gaussian_cloud1.numPoints))
    max_rot_diff_idx = np.argmax(rot_diff)
    test_case.assertTrue(
        np.allclose(rotations_dot, np.ones(gaussian_cloud1.numPoints), atol=MAX_ROTATION_ERROR),
        f"Max rotation diff {rot_diff[max_rot_diff_idx]} at index {max_rot_diff_idx}: quaternion {q1s[max_rot_diff_idx]} vs {q2s[max_rot_diff_idx]}",
    )

    alphas1 = sigmoid(np.array(gaussian_cloud1.alphas))
    alphas2 = sigmoid(np.array(gaussian_cloud2.alphas))
    alpha_diff = np.abs(alphas1 - alphas2)
    max_alpha_diff_idx = np.argmax(alpha_diff)
    test_case.assertTrue(
        np.allclose(alphas1, alphas2, atol=MAX_ALPHA_ERROR),
        f"Max alpha diff {alpha_diff[max_alpha_diff_idx]} at index {max_alpha_diff_idx}: {alphas1[max_alpha_diff_idx]} vs {alphas2[max_alpha_diff_idx]}",
    )

    gaussian_cloud1_colors = (np.array(gaussian_cloud1.colors).reshape(-1, 3) * COLOR_SCALE * 255.0 + (0.5 * 255.0)).clip(0, 255)
    gaussian_cloud2_colors = (np.array(gaussian_cloud2.colors).reshape(-1, 3) * COLOR_SCALE * 255.0 + (0.5 * 255.0)).clip(0, 255)
    color_diff = np.abs(gaussian_cloud1_colors - gaussian_cloud2_colors)
    max_color_diff_idx = np.unravel_index(np.argmax(color_diff), color_diff.shape)
    test_case.assertTrue(
        np.allclose(gaussian_cloud1_colors, gaussian_cloud2_colors, atol=MAX_COLOR_ERROR),
        f"Max color diff {color_diff[max_color_diff_idx]} at index {max_color_diff_idx}: RGB {gaussian_cloud1_colors[max_color_diff_idx[0]]} vs {gaussian_cloud2_colors[max_color_diff_idx[0]]}",
    )
    if len(gaussian_cloud1.sh) > 0:
        sh_diff = np.abs(np.array(gaussian_cloud1.sh) - np.array(gaussian_cloud2.sh))
        max_sh_diff_idx = np.argmax(sh_diff)
        test_case.assertTrue(
            np.allclose(gaussian_cloud1.sh, gaussian_cloud2.sh, atol=MAX_SH_ERROR),
            f"Max SH diff {sh_diff[max_sh_diff_idx]} at index {max_sh_diff_idx}: {gaussian_cloud1.sh[max_sh_diff_idx]} vs {gaussian_cloud2.sh[max_sh_diff_idx]}",
        )


class TestPly(SpzTestCase):
    def test_ply2spz_flywheel_float32_sh3(self):
        gt_file = self.input_test_path() / "ply/shoe_sh3_float32.ply"
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(gt_file)
        spz_file = self.output_test_path() / "shoe_sh3_float32.spz"
        gaussian_cloud_to_spz_file(gaussian_cloud_ply, spz_file, disable_sh_min_max_scaling=True)
        gaussian_cloud_spz = read_spz_to_gaussian_cloud(spz_file)
        gaussian_cloud_to_ply_file(gaussian_cloud_spz, self.output_test_path() / "shoe_sh3_float32_reconstructed.ply")
        _assert_gaussian_cloud_equal(self, gaussian_cloud_ply, gaussian_cloud_spz)

    def test_ply2spz_flywheel_uint16_sh0(self):
        gt_file = self.input_test_path() / "ply/shoe_sh0_uint16.ply"
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(gt_file)
        spz_file = self.output_test_path() / "shoe_sh0_uint16.spz"
        gaussian_cloud_to_spz_file(gaussian_cloud_ply, spz_file)
        gaussian_cloud_spz = read_spz_to_gaussian_cloud(spz_file)
        gaussian_cloud_to_ply_file(gaussian_cloud_spz, self.output_test_path() / "shoe_sh0_uint16_reconstructed.ply")
        _assert_gaussian_cloud_equal(self, gaussian_cloud_ply, gaussian_cloud_spz)

    def test_ply2spz_flywheel_uint16_sh3(self):
        gt_file = self.input_test_path() / "ply/tree_sh3_uint16.ply"
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(gt_file)
        spz_file = self.output_test_path() / "tree_sh3_uint16.spz"
        gaussian_cloud_to_spz_file(gaussian_cloud_ply, spz_file, disable_sh_min_max_scaling=True)
        gaussian_cloud_spz = read_spz_to_gaussian_cloud(spz_file)
        gaussian_cloud_to_ply_file(gaussian_cloud_spz, self.output_test_path() / "tree_sh3_uint16_reconstructed.ply")
        _assert_gaussian_cloud_equal(self, gaussian_cloud_ply, gaussian_cloud_spz)

    def test_ply2spz_float16_sh3_toy(self):
        gt_file = self.input_test_path() / "ply/Furry_Toy_2_sh3_float16.ply"
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(gt_file)
        spz_file = self.output_test_path() / "Furry_Toy_2_sh3_float16.spz"
        gaussian_cloud_to_spz_file(gaussian_cloud_ply, spz_file)
        gaussian_cloud_spz = read_spz_to_gaussian_cloud(spz_file)
        gaussian_cloud_to_ply_file(gaussian_cloud_spz, self.output_test_path() / "Furry_Toy_2_sh3_float16_reconstructed.ply")
        _assert_gaussian_cloud_equal(self, gaussian_cloud_ply, gaussian_cloud_spz)

    def test_ply2spz_float32_sh3_toy(self):
        gt_file = self.input_test_path() / "ply/Furry_Toy_2_sh3_float32.ply"
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(gt_file)
        spz_file = self.output_test_path() / "Furry_Toy_2_sh3_float32.spz"
        gaussian_cloud_to_spz_file(gaussian_cloud_ply, spz_file)
        gaussian_cloud_spz = read_spz_to_gaussian_cloud(spz_file)
        gaussian_cloud_to_ply_file(gaussian_cloud_spz, self.output_test_path() / "Furry_Toy_2_sh3_float32_reconstructed.ply")
        _assert_gaussian_cloud_equal(self, gaussian_cloud_ply, gaussian_cloud_spz)

    def test_ply2spz_neuralcg_float32_sh4(self):
        gt_file = self.input_test_path() / "ply/plant_sh4_float32.ply"
        gaussian_cloud_ply = read_ply_to_gaussian_cloud(gt_file)
        spz_file = self.output_test_path() / "plant_sh4_float32.spz"
        gaussian_cloud_to_spz_file(gaussian_cloud_ply, spz_file)
        gaussian_cloud_spz = read_spz_to_gaussian_cloud(spz_file)
        gaussian_cloud_to_ply_file(gaussian_cloud_spz, self.output_test_path() / "plant_sh4_float32_reconstructed.ply")
        _assert_gaussian_cloud_equal(self, gaussian_cloud_ply, gaussian_cloud_spz)
