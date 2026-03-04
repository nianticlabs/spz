"""Pytest conftest.py - shared fixtures and pytest configuration."""

import sys
from pathlib import Path

# Add the tests/python directory to the Python path so test_utils can be imported
sys.path.insert(0, str(Path(__file__).parent))

