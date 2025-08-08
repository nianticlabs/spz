from spz.spz_test_case import SpzTestCase


class TestDataDownload(SpzTestCase):

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
