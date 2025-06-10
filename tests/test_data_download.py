# This file is adapted from https://git.corp.adobe.com/3di/python-scaffold

from spz.spz_test_case import SpzTestCase


class TestDataDownload(SpzTestCase):
    def test_downloaded_npz(self):
        ut_file = self.input_test_path() / "npz/v0.4/Furry_Toy_2.npz"
        sha_file = self.file_sha256(ut_file)
        self.assertEqual(sha_file, "1ec24f805130fef5807f8f91f5ae7fb7dd699276ce3a853a9c3746ee7c411d2f")

    def test_downloaded_gltf(self):
        ut_file = self.input_test_path() / "gltf/v0.4/Furry_Toy_2.gltf"
        sha_file = self.file_sha256(ut_file)
        self.assertEqual(sha_file, "305c5b7a33c4be31e33f1c554ae6e27cf533a6f7c1636d767cec422d34a5454d")

    def test_downloaded_usdc(self):
        ut_file = self.input_test_path() / "usdc/v0.4/Furry_Toy_2.usdc"
        sha_file = self.file_sha256(ut_file)
        self.assertEqual(sha_file, "95d7319f9d05cfc1f7eaa0f3c9f60283b259674763fd3fe9f582ca092592049f")

    def test_downloaded_gs_f16(self):
        ut_file = self.input_test_path() / "ply/Furry_Toy_2_sh3_float16.ply"
        sha_file = self.file_sha256(ut_file)
        self.assertEqual(sha_file, "41dfb5d9c2ad224535a784a87f9923b9eda6aedacb436685245b2fe4c17c5986")

    def test_downloaded_gs_f32(self):
        ut_file = self.input_test_path() / "ply/Furry_Toy_2_sh3_float32.ply"
        sha_file = self.file_sha256(ut_file)
        self.assertEqual(sha_file, "807aae566a75f4b6c00d6c83c430805421ced004979812dcd5bd38ebf50ec688")

    def test_downloaded_gs_sh3_f16(self):
        ut_file = self.input_test_path() / "ply/shoe_sh3_float32.ply"
        sha_file = self.file_sha256(ut_file)
        self.assertEqual(sha_file, "2d47198c957630638999483a2bcc91205815d16e27510a2f2fb50ccb4d4b638d")

    def test_downloaded_gs_sh0_uint16(self):
        ut_file = self.input_test_path() / "ply/shoe_sh0_uint16.ply"
        sha_file = self.file_sha256(ut_file)
        self.assertEqual(sha_file, "f61d99c8b63205da4e5670cac12b901a808804ec9b1f1fe1dc574d59edea6a8b")

    def test_downloaded_gs_sh3_uint16(self):
        ut_file = self.input_test_path() / "ply/tree_sh3_uint16.ply"
        sha_file = self.file_sha256(ut_file)
        self.assertEqual(sha_file, "3d453904e4abc7bd4af9dac4d94a6dd85321a1263bee626d1d2a80c36a4feee3")
