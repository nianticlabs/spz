"""Setup script for spz with automatic C++ bindings build.

This setup script automatically builds the SPZ Python bindings from C++ code
using Metabuild during pip install.
"""

import platform
import subprocess
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class MetabuildExtension(Extension):
    """Metabuild Extension to specify the source directory."""

    def __init__(self, name, sourcedir="."):
        super().__init__(name, sources=[])
        self.sourcedir = Path(sourcedir).resolve()


class MetabuildSpzBindings(build_ext):
    """Metabuild build_ext command for building SPZ bindings."""

    def run(self):
        # Ensure Metabuild is installed
        try:
            subprocess.check_call(["metabuild", "version"])
        except OSError as exc:
            raise RuntimeError("Metabuild must be installed to build this extension") from exc

        # Build spz library
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        build_temp = Path(self.build_temp).resolve()

        # Prepare and build using Metabuild
        build_temp.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(["metabuild", "prepare", "--verbose"], cwd=ext.sourcedir)
        subprocess.check_call(["metabuild", "build", "-c", "Release", "--verbose"], cwd=ext.sourcedir)

        # Determine the output file based on the platform
        system = platform.system()
        if system == "Darwin":  # macOS
            output_file = ext.sourcedir / "build_mb/xcode_macos/Release/universal/build/spz/spz/spz_bindings_so/spz_bindings.so"
        elif system == "Windows":
            output_file = ext.sourcedir / "build_mb/msvs_win32/Release/x64/build/spz/spz/spz_bindings_pyd/spz_bindings.pyd"
        elif system == "Linux":
            output_file = ext.sourcedir / "build_mb/cmake_linux/Release/x64/build/spz/spz/spz_bindings_so/spz_bindings.so"
        else:
            raise RuntimeError(f"Unsupported platform: {system}")

        # Check if the output file exists
        if not output_file.exists():
            raise RuntimeError(f"Expected shared library not found: {output_file}")

        # Copy the resulting file to the local spz directory
        ext_filename = Path(self.get_ext_fullpath(ext.name)).name
        target_dir = Path(ext.sourcedir) / "spz"
        target_file = target_dir / ext_filename

        # Ensure the target directory exists
        target_dir.mkdir(parents=True, exist_ok=True)
        self.copy_file(str(output_file), str(target_file))


if __name__ == "__main__":
    setup(
        ext_modules=[MetabuildExtension("spz.spz_bindings", sourcedir=".")],
        cmdclass={"build_ext": MetabuildSpzBindings},
        include_package_data=True,
        package_data={
            "spz.spz": ["*.so", "*.pyd"],
        },
    )
