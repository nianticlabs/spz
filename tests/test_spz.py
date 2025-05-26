import pytest
import os
import tempfile
import struct
import spz


def test_import():
    """Test that the module can be imported and has expected attributes."""
    assert hasattr(spz, 'ply_to_spz')
    assert hasattr(spz, '__version__')


def test_ply_to_spz_missing_file():
    """Test that conversion fails gracefully with missing input file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ply_path = os.path.join(tmpdir, "nonexistent.ply")
        spz_path = os.path.join(tmpdir, "output.spz")
        
        result = spz.ply_to_spz(ply_path, spz_path)
        assert result is False


def test_coordinate_systems_nonexistent_file():
    """Test that different coordinate systems are accepted without error."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ply_path = os.path.join(tmpdir, "nonexistent.ply")
        spz_path = os.path.join(tmpdir, "output.spz")
        
        # These should all fail gracefully (return False) but not raise exceptions
        coordinate_systems = ["RDF", "RUB", "LUF", "RUF", "UNSPECIFIED"]
        
        for cs in coordinate_systems:
            result = spz.ply_to_spz(ply_path, spz_path, coordinate_system=cs)
            assert result is False  # Should fail because file doesn't exist


def test_invalid_coordinate_system():
    """Test that invalid coordinate systems raise appropriate errors."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ply_path = os.path.join(tmpdir, "nonexistent.ply")
        spz_path = os.path.join(tmpdir, "output.spz")
        
        with pytest.raises(Exception):  # Should raise an exception for invalid CS
            spz.ply_to_spz(ply_path, spz_path, coordinate_system="INVALID")


def create_simple_ply(ply_path):
    """Create a simple PLY file with correct format for testing."""
    # Create a minimal PLY file with 2 gaussians in the expected format
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
property float f_rest_0
property float f_rest_1
property float f_rest_2
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
    
    # Create some simple gaussian data (2 gaussians with 18 floats each)
    gaussian1 = [
        0.0, 0.0, 0.0,  # position (x, y, z)
        0.0, 0.0, 1.0,  # normal (nx, ny, nz) - not really used
        0.5, 0.3, 0.2,  # color DC (f_dc_0, f_dc_1, f_dc_2)
        0.0, 0.0, 0.0,  # SH coeffs (f_rest_0, f_rest_1, f_rest_2)
        1.0,            # opacity
        0.1, 0.1, 0.1,  # scale (scale_0, scale_1, scale_2)
        1.0, 0.0, 0.0, 0.0  # rotation quaternion (rot_0, rot_1, rot_2, rot_3)
    ]
    
    gaussian2 = [
        1.0, 1.0, 1.0,  # position (x, y, z)
        0.0, 0.0, 1.0,  # normal (nx, ny, nz)
        0.8, 0.1, 0.1,  # color DC (f_dc_0, f_dc_1, f_dc_2)
        0.0, 0.0, 0.0,  # SH coeffs (f_rest_0, f_rest_1, f_rest_2)
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


def test_ply_to_spz_simple_conversion():
    """Test PLY to SPZ conversion with a simple compatible PLY file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ply_path = os.path.join(tmpdir, "simple.ply")
        spz_path = os.path.join(tmpdir, "simple.spz")
        
        # Create a simple PLY file in the correct format
        create_simple_ply(ply_path)
        
        # Test conversion
        result = spz.ply_to_spz(ply_path, spz_path, coordinate_system="RDF")
        
        # Verify conversion succeeded
        assert result is True, "Simple PLY to SPZ conversion should succeed"
        
        # Verify output file was created
        assert os.path.exists(spz_path), "SPZ output file should be created"
        
        # Verify output file has content
        spz_size = os.path.getsize(spz_path)
        assert spz_size > 0, "SPZ file should not be empty"
        
        # Verify compression
        ply_size = os.path.getsize(ply_path)
        
        print(f"Simple conversion test passed:")
        print(f"  PLY size: {ply_size:,} bytes")
        print(f"  SPZ size: {spz_size:,} bytes")


def test_ply_to_spz_actual_conversion():
    """Test actual PLY to SPZ conversion using the sample file."""
    # Get the path to the sample PLY file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    ply_path = os.path.join(project_root, "samples", "splat.ply")
    
    # Skip test if sample file doesn't exist
    if not os.path.exists(ply_path):
        pytest.skip("Sample PLY file not found")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        spz_path = os.path.join(tmpdir, "output.spz")
        
        # Test conversion with default coordinate system (RDF)
        result = spz.ply_to_spz(ply_path, spz_path, coordinate_system="RDF")
        
        # The PLY file should now have the correct float format and conversion should succeed
        assert result is True, "PLY to SPZ conversion should succeed with corrected float format"
        
        # Verify output file was created
        assert os.path.exists(spz_path), "SPZ output file should be created"
        
        # Verify output file has content
        spz_size = os.path.getsize(spz_path)
        assert spz_size > 0, "SPZ file should not be empty"
        
        # Verify compression
        ply_size = os.path.getsize(ply_path)
        compression_ratio = ply_size / spz_size
        
        print(f"Actual PLY conversion test passed:")
        print(f"  PLY size: {ply_size:,} bytes")
        print(f"  SPZ size: {spz_size:,} bytes") 
        print(f"  Compression ratio: {compression_ratio:.2f}x")
        
        # SPZ should be significantly smaller than PLY
        assert compression_ratio > 1.0, "SPZ file should be smaller than PLY file"


def test_ply_to_spz_different_coordinate_systems():
    """Test PLY to SPZ conversion with different coordinate systems."""
    # Get the path to the sample PLY file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    ply_path = os.path.join(project_root, "samples", "splat.ply")
    
    # Skip test if sample file doesn't exist
    if not os.path.exists(ply_path):
        pytest.skip("Sample PLY file not found")
    
    coordinate_systems = ["RDF", "RUB", "LUF", "RUF"]
    
    with tempfile.TemporaryDirectory() as tmpdir:
        for cs in coordinate_systems:
            spz_path = os.path.join(tmpdir, f"output_{cs}.spz")
            
            # Test conversion with each coordinate system
            result = spz.ply_to_spz(ply_path, spz_path, coordinate_system=cs)
            
            # With the corrected PLY format, conversion should succeed for all coordinate systems
            assert result is True, f"PLY to SPZ conversion should succeed with {cs} coordinate system"
            
            # Verify the output
            assert os.path.exists(spz_path), f"SPZ file with {cs} coordinate system should be created"
            assert os.path.getsize(spz_path) > 0, f"SPZ file with {cs} coordinate system should not be empty"
            print(f"Conversion with {cs} coordinate system succeeded") 