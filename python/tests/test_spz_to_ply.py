"""Tests for SPZ to PLY conversion functionality."""

import pytest
import os
import tempfile
import spz
from .conftest import create_simple_ply


def test_spz_to_ply_basic(sample_ply_path):
    """Test basic SPZ to PLY conversion functionality."""
    # Skip test if sample file doesn't exist
    if not os.path.exists(sample_ply_path):
        pytest.skip("Sample PLY file not found")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        spz_path = os.path.join(tmpdir, "intermediate.spz")
        output_ply_path = os.path.join(tmpdir, "output.ply")
        
        # First convert PLY to SPZ
        result1 = spz.ply_to_spz(sample_ply_path, spz_path, coordinate_system="RDF")
        assert result1, "PLY to SPZ conversion should succeed"
        assert os.path.exists(spz_path), "SPZ file should be created"
        
        # Then convert SPZ back to PLY
        result2 = spz.spz_to_ply(spz_path, output_ply_path, coordinate_system="RDF")
        assert result2, "SPZ to PLY conversion should succeed"
        assert os.path.exists(output_ply_path), "Output PLY file should be created"
        assert os.path.getsize(output_ply_path) > 0, "Output PLY file should not be empty"
        
        # Check that we achieved compression in the intermediate SPZ file
        original_size = os.path.getsize(sample_ply_path)
        spz_size = os.path.getsize(spz_path)
        output_size = os.path.getsize(output_ply_path)
        
        compression_ratio = original_size / spz_size
        print(f"Original PLY: {original_size:,} bytes")
        print(f"SPZ: {spz_size:,} bytes")
        print(f"Output PLY: {output_size:,} bytes")
        print(f"Compression ratio (original/SPZ): {compression_ratio:.2f}x")
        
        # SPZ should be significantly smaller than PLY
        assert compression_ratio > 5, "SPZ should provide significant compression"


def test_spz_to_ply_different_coordinate_systems(sample_ply_path):
    """Test SPZ to PLY conversion with different coordinate systems."""
    # Skip test if sample file doesn't exist
    if not os.path.exists(sample_ply_path):
        pytest.skip("Sample PLY file not found")
    
    coordinate_systems = ["RDF", "RUB", "LUF", "RUF", "UNSPECIFIED"]
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # First create an SPZ file to work with
        spz_path = os.path.join(tmpdir, "test.spz")
        result = spz.ply_to_spz(sample_ply_path, spz_path, coordinate_system="RDF")
        assert result, "PLY to SPZ conversion should succeed"
        
        # Test conversion to PLY with each coordinate system
        for cs in coordinate_systems:
            output_ply_path = os.path.join(tmpdir, f"output_{cs}.ply")
            
            result = spz.spz_to_ply(spz_path, output_ply_path, coordinate_system=cs)
            assert result, f"SPZ to PLY conversion with {cs} coordinate system should succeed"
            assert os.path.exists(output_ply_path), f"PLY file with {cs} coordinate system should be created"
            assert os.path.getsize(output_ply_path) > 0, f"PLY file with {cs} coordinate system should not be empty"
            print(f"SPZ to PLY conversion with {cs} succeeded")


def test_spz_to_ply_error_handling(sample_ply_path):
    """Test error handling for spz_to_ply function."""
    with tempfile.TemporaryDirectory() as tmpdir:
        nonexistent_spz = os.path.join(tmpdir, "nonexistent.spz")
        output_ply = os.path.join(tmpdir, "output.ply")
        
        # Test with nonexistent SPZ file
        result = spz.spz_to_ply(nonexistent_spz, output_ply)
        assert not result, "Should return False for nonexistent SPZ file"
        
        # Test with invalid coordinate system
        if os.path.exists(sample_ply_path):
            valid_spz = os.path.join(tmpdir, "valid.spz")
            spz.ply_to_spz(sample_ply_path, valid_spz)
            
            with pytest.raises(ValueError, match="Invalid coordinate system"):
                spz.spz_to_ply(valid_spz, output_ply, coordinate_system="INVALID")


def test_spz_to_ply_simple_conversion():
    """Test SPZ to PLY conversion with a simple created PLY file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ply_path = os.path.join(tmpdir, "simple.ply")
        spz_path = os.path.join(tmpdir, "simple.spz")
        output_ply_path = os.path.join(tmpdir, "output.ply")
        
        # Create a simple PLY file
        create_simple_ply(ply_path)
        
        # Convert PLY to SPZ
        result1 = spz.ply_to_spz(ply_path, spz_path, coordinate_system="RDF")
        assert result1, "PLY to SPZ conversion should succeed"
        
        # Convert SPZ back to PLY
        result2 = spz.spz_to_ply(spz_path, output_ply_path, coordinate_system="RDF")
        assert result2, "SPZ to PLY conversion should succeed"
        
        # Check files exist and have content
        assert os.path.exists(output_ply_path), "Output PLY file should be created"
        assert os.path.getsize(output_ply_path) > 0, "Output PLY file should not be empty"
        
        # Check compression was achieved
        ply_size = os.path.getsize(ply_path)
        spz_size = os.path.getsize(spz_path)
        output_size = os.path.getsize(output_ply_path)
        
        print(f"Simple SPZ to PLY conversion test:")
        print(f"  Original PLY: {ply_size:,} bytes")
        print(f"  SPZ: {spz_size:,} bytes")
        print(f"  Output PLY: {output_size:,} bytes")
        print(f"  Compression ratio: {ply_size/spz_size:.2f}x")
        
        assert spz_size < ply_size, "SPZ should be smaller than original PLY" 