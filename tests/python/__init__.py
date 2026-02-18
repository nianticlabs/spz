# Test package for spz Python bindings

# On Windows, configure DLL search paths before importing spz
# This is required for stable ABI wheels that link against python3.dll
import os
import sys
from pathlib import Path

if sys.platform == "win32":
    # Add the base Python installation directory (where python3.dll lives)
    python_base = Path(sys.base_prefix)
    if python_base.exists():
        os.add_dll_directory(str(python_base))

    # Add the spz package directory to find any built DLLs (for wheel installs)
    spz_package_dir = Path(__file__).parent.parent.parent / "src" / "python" / "spz"
    if spz_package_dir.exists():
        os.add_dll_directory(str(spz_package_dir))

    # For editable installs, also add the build directory where zlib.dll is located
    # The build directory pattern is build/{wheel_tag}/_deps/zlib-build/Release
    build_dir = Path(__file__).parent.parent.parent / "build"
    if build_dir.exists():
        # Find zlib.dll in the build directory
        for zlib_dll_dir in build_dir.glob("*/_deps/zlib-build/Release"):
            if zlib_dll_dir.exists():
                os.add_dll_directory(str(zlib_dll_dir.resolve()))

# Export shared test utilities
from .test_utils import (
    sh_epsilon,
    normalized,
    axis_angle_quat,
    times,
    read_file,
    make_test_gaussian_cloud,
)

