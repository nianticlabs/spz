"""
SPZ: Python bindings for compressed 3D Gaussian splats

This package provides Python bindings for loading, converting, and saving
SPZ files, which are a compressed format for 3D Gaussian splats.
"""

from .core import ply_to_spz

__version__ = "1.1.0"
__all__ = ["ply_to_spz"] 