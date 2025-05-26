"""Shared test fixtures and utilities for SPZ tests."""

import pytest
import os
import struct
import tempfile


@pytest.fixture
def sample_ply_path():
    """Get the path to the sample PLY file."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(current_dir))
    return os.path.join(project_root, "samples", "splat.ply")


def create_simple_ply(ply_path):
    """Create a simple PLY file with correct format for testing (degree 0 SH - no extra coefficients)."""
    # Create a minimal PLY file with 2 gaussians in degree 0 format (no spherical harmonics beyond DC)
    header = """ply
format binary_little_endian 1.0
element vertex 2
property float x
property float y
property float z
property float nx
property float ny
property float nz
property float f_dc_0
property float f_dc_1
property float f_dc_2
property float opacity
property float scale_0
property float scale_1
property float scale_2
property float rot_0
property float rot_1
property float rot_2
property float rot_3
end_header
"""
    
    # Create some simple gaussian data (2 gaussians with degree 0 SH - only DC coefficients)
    gaussian1 = [
        0.0, 0.0, 0.0,  # position (x, y, z)
        0.0, 0.0, 1.0,  # normal (nx, ny, nz) - not really used
        0.5, 0.3, 0.2,  # color DC (f_dc_0, f_dc_1, f_dc_2)
        1.0,            # opacity
        0.1, 0.1, 0.1,  # scale (scale_0, scale_1, scale_2)
        1.0, 0.0, 0.0, 0.0  # rotation quaternion (rot_0, rot_1, rot_2, rot_3)
    ]
    
    gaussian2 = [
        1.0, 1.0, 1.0,  # position (x, y, z)
        0.0, 0.0, 1.0,  # normal (nx, ny, nz)
        0.8, 0.1, 0.1,  # color DC (f_dc_0, f_dc_1, f_dc_2)
        0.8,            # opacity
        0.2, 0.15, 0.1, # scale (scale_0, scale_1, scale_2)
        0.707, 0.707, 0.0, 0.0  # rotation quaternion (rot_0, rot_1, rot_2, rot_3)
    ]
    
    with open(ply_path, 'wb') as f:
        f.write(header.encode('ascii'))
        
        # Write binary data for both gaussians
        for gaussian in [gaussian1, gaussian2]:
            for value in gaussian:
                f.write(struct.pack('<f', value))  # little-endian float 