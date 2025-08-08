"""Setup script for spz with automatic C++ bindings build.

This setup script automatically builds the SPZ Python bindings from C++ code
using Metabuild during pip install.
"""

import platform
import subprocess
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    """CMake Extension to specify the source directory."""

    def __init__(self, name, sourcedir="."):
        super().__init__(name, sources=[])
        self.sourcedir = Path(sourcedir).resolve()


class CMakeSpzBindings(build_ext):
    """CMake build_ext command for building SPZ bindings."""

    def run(self):
        # Ensure CMake is installed
        try:
            subprocess.check_call(["cmake", "--version"])
        except OSError as exc:
            raise RuntimeError("CMake must be installed to build this extension") from exc

        # Build spz library
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        build_temp = Path(self.build_temp).resolve()
        build_temp.mkdir(parents=True, exist_ok=True)
        config = "Release"

        # Run Metabuild steps
        subprocess.check_call([
            "cmake",
            "-B",
            "build",
            "-S",
            ".",
            "-DSPZ_BUILD_PYTHON_BINDINGS=ON",
            f"-DCMAKE_BUILD_TYPE={config}"
        ], cwd=ext.sourcedir)
        subprocess.check_call(["cmake", "--build", "build", "--config", config], cwd=ext.sourcedir)

        # Determine the output file
        system = platform.system()
        if system == "Darwin":
            output_file = ext.sourcedir / "build" / "spz_bindings.so"
        elif system == "Windows":
            output_file = ext.sourcedir / "build" / config / "spz_bindings.pyd"
        elif system == "Linux":
            output_file = ext.sourcedir / "build" / "spz_bindings.so"
        else:
            raise RuntimeError(f"Unsupported platform: {system}")

        if not output_file.exists():
            raise RuntimeError(f"Expected shared library not found: {output_file}")

        # Copy into setuptools' expected build directory
        target_path = Path(self.get_ext_fullpath(ext.name))
        target_path.parent.mkdir(parents=True, exist_ok=True)
        self.copy_file(str(output_file), str(target_path))


if __name__ == "__main__":
    setup(
        ext_modules=[CMakeExtension("spz.spz_bindings", sourcedir=".")],
        cmdclass={"build_ext": CMakeSpzBindings},
        include_package_data=True,
        package_data={
            "spz.spz": ["*.so", "*.pyd"],
        },
    )
