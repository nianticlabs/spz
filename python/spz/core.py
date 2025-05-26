"""
Core Python interface for SPZ operations.
"""

from typing import Optional
from ._core import _ply_to_spz_impl


def ply_to_spz(
    ply_path: str,
    spz_path: str,
    coordinate_system: str = "RDF"
) -> bool:
    """
    Convert a PLY file to SPZ format.
    
    Args:
        ply_path: Path to the input PLY file
        spz_path: Path to the output SPZ file
        coordinate_system: Coordinate system of the input PLY file.
                          Options: "RDF" (default, typical PLY format),
                                  "RUB" (Three.js/OpenGL),
                                  "LUF" (GLB format), 
                                  "RUF" (Unity format),
                                  or "UNSPECIFIED"
    
    Returns:
        True if conversion was successful, False otherwise
        
    Example:
        >>> import spz
        >>> success = spz.ply_to_spz("input.ply", "output.spz")
        >>> if success:
        ...     print("Conversion successful!")
    """
    return _ply_to_spz_impl(ply_path, spz_path, coordinate_system) 