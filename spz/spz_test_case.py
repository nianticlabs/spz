import hashlib
import logging
import unittest
from pathlib import Path

from spz.log import Channel, logger


class SpzTestCase(unittest.TestCase):
    root_dir = Path(__file__).absolute().parents[1]

    @classmethod
    def setUpClass(cls):
        logger.info("setUpClass %r", cls.__name__)
        if len(logger.handlers) == 0:
            logger.setLevel(Channel.INFO)
            print_to_screen = logging.StreamHandler()
            print_to_screen.setLevel(Channel.INFO)
            logger.addHandler(print_to_screen)

    @classmethod
    def tearDownClass(cls):
        logger.info("tearDown %r", cls.__name__)

    @staticmethod
    def sha256(some_value):
        if isinstance(some_value, str):
            some_value = some_value.encode("utf-8")
        hash_sha256 = hashlib.sha256()
        hash_sha256.update(some_value)
        return hash_sha256.hexdigest()

    def read_file(self, path: Path, asbytes=False):
        self.assertTrue(path.is_file())
        if asbytes:
            return path.read_bytes()
        return path.read_text(encoding="UTF-8")

    def file_sha256(self, path: Path) -> str:
        return self.sha256(self.read_file(path=path, asbytes=True))

    @classmethod
    def input_test_path(cls):
        return cls.root_dir / "data"

    @classmethod
    def output_test_path(cls):
        return cls.root_dir / ".tmp_test_out"
