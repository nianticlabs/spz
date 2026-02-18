"""Pytest conftest.py - shared fixtures and pytest configuration."""

import os
import sys
from pathlib import Path

# On Windows, add Python base directory to DLL search path for stable ABI support
# This MUST be done before any imports of the spz module
# This ensures python3.dll can be found when loading the spz extension module
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

# Add the tests/python directory to the Python path so test_utils can be imported
sys.path.insert(0, str(Path(__file__).parent))

