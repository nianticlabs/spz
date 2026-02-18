# Test package for spz Python bindings

# Configure Windows DLL search paths before importing spz
from pathlib import Path
from ._dll_setup import configure_dll_search_paths
configure_dll_search_paths(Path(__file__).parent.parent.parent)

# Export shared test utilities
from .test_utils import (
    sh_epsilon,
    normalized,
    axis_angle_quat,
    times,
    read_file,
    make_test_gaussian_cloud,
)

