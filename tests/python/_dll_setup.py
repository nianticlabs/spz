"""Configure Windows DLL search paths for stable ABI support.

This must be called before any imports of the spz module.
Stable ABI wheels (cp312 abi3) link against python3.dll, which may not be
on the default DLL search path. This helper adds the necessary directories.
"""

import os
import sys
from pathlib import Path


def configure_dll_search_paths(root_dir):
    """Add DLL search directories needed to load spz on Windows.

    Args:
        root_dir: The project root directory (typically 3 levels up from the
                  calling file in tests/python/).
    """
    if sys.platform != "win32":
        return

    root_dir = Path(root_dir)

    # Add the base Python installation directory (where python3.dll lives)
    python_base = Path(sys.base_prefix)
    if python_base.exists():
        os.add_dll_directory(str(python_base))

    # Add the spz package directory to find any built DLLs (for wheel installs)
    spz_package_dir = root_dir / "src" / "python" / "spz"
    if spz_package_dir.exists():
        os.add_dll_directory(str(spz_package_dir))

    # For editable installs, also add the build directory where zlib.dll is located
    # The build directory pattern is build/{wheel_tag}/_deps/zlib-build/Release
    build_dir = root_dir / "build"
    if build_dir.exists():
        for zlib_dll_dir in build_dir.glob("*/_deps/zlib-build/Release"):
            if zlib_dll_dir.exists():
                os.add_dll_directory(str(zlib_dll_dir.resolve()))
