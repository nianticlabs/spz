"""
Core Python interface for SPZ operations.
"""

from typing import Optional
from ._core import _ply_to_spz_impl, _spz_to_ply_impl


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


def spz_to_ply(
    spz_path: str,
    ply_path: str,
    coordinate_system: str = "RDF"
) -> bool:
    """
    Convert an SPZ file to PLY format.
    
    Args:
        spz_path: Path to the input SPZ file
        ply_path: Path to the output PLY file  
        coordinate_system: Target coordinate system for the output PLY file.
                          Options: "RDF" (default, typical PLY format),
                                  "RUB" (Three.js/OpenGL),
                                  "LUF" (GLB format),
                                  "RUF" (Unity format),
                                  or "UNSPECIFIED"
    
    Returns:
        True if conversion was successful, False otherwise
        
    Example:
        >>> import spz
        >>> success = spz.spz_to_ply("input.spz", "output.ply")
        >>> if success:
        ...     print("Conversion successful!")
    """
    return _spz_to_ply_impl(spz_path, ply_path, coordinate_system) 