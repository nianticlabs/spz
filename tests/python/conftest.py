"""Pytest conftest.py - shared fixtures and pytest configuration."""

import importlib.util
import sys
from pathlib import Path

# Configure Windows DLL search paths before any imports of the spz module.
# We load _dll_setup directly by path since conftest.py is not part of a
# package and tests/python/ is not yet on sys.path at this point.
_dll_setup_path = Path(__file__).parent / "_dll_setup.py"
_spec = importlib.util.spec_from_file_location("_dll_setup", _dll_setup_path)
_dll_setup = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_dll_setup)
_dll_setup.configure_dll_search_paths(Path(__file__).parent.parent.parent)

# Add the tests/python directory to the Python path so test_utils can be imported
sys.path.insert(0, str(Path(__file__).parent))

