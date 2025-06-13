# This file is adapted from https://git.corp.adobe.com/euclid/python-project-scaffold
# This is for storing your large data outside of git.

import pathlib
import sys

from tiny.commands.main import main as tiny_main


def main(yamlpath=sys.argv[1], operation=sys.argv[2]):
    config_file = pathlib.Path(yamlpath)
    if not config_file.is_file():
        config_file = pathlib.Path(__file__).parent / "tiny.yaml"
    if not config_file.is_file():
        raise BaseException("Expected an existing YAML manifest file as the first argument")

    arg_array = [
        operation,
        "--num-retry",
        "10",  # try 10 times.
        "--retry-wait",
        "1",  # wait for 1 second after failure
        "--num-sim-jobs",
        "20",  # how many files to process at once
        "--config",
        str(config_file),
    ]

    print(f"Running `tiny {' '.join(arg_array)}` ...")
    tiny_main(arg_array)


if __name__ == "__main__":
    main()
