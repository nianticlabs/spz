"""Tests for PLY to SPZ conversion functionality."""

import pytest
import os
import tempfile
import spz
from .conftest import create_simple_ply


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


def test_ply_to_spz_actual_conversion(sample_ply_path):
    """Test actual PLY to SPZ conversion using the sample file."""
    # Skip test if sample file doesn't exist
    if not os.path.exists(sample_ply_path):
        pytest.skip("Sample PLY file not found")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        spz_path = os.path.join(tmpdir, "output.spz")
        
        # Test conversion with default coordinate system (RDF)
        result = spz.ply_to_spz(sample_ply_path, spz_path, coordinate_system="RDF")
        
        # The PLY file should now have the correct float format and conversion should succeed
        assert result is True, "PLY to SPZ conversion should succeed with corrected float format"
        
        # Verify output file was created
        assert os.path.exists(spz_path), "SPZ output file should be created"
        
        # Verify output file has content
        spz_size = os.path.getsize(spz_path)
        assert spz_size > 0, "SPZ file should not be empty"
        
        # Verify compression
        ply_size = os.path.getsize(sample_ply_path)
        compression_ratio = ply_size / spz_size
        
        print(f"Actual PLY conversion test passed:")
        print(f"  PLY size: {ply_size:,} bytes")
        print(f"  SPZ size: {spz_size:,} bytes") 
        print(f"  Compression ratio: {compression_ratio:.2f}x")
        
        # SPZ should be significantly smaller than PLY
        assert compression_ratio > 1.0, "SPZ file should be smaller than PLY file"


def test_ply_to_spz_different_coordinate_systems(sample_ply_path):
    """Test PLY to SPZ conversion with different coordinate systems."""
    # Skip test if sample file doesn't exist
    if not os.path.exists(sample_ply_path):
        pytest.skip("Sample PLY file not found")
    
    coordinate_systems = ["RDF", "RUB", "LUF", "RUF"]
    
    with tempfile.TemporaryDirectory() as tmpdir:
        for cs in coordinate_systems:
            spz_path = os.path.join(tmpdir, f"output_{cs}.spz")
            
            # Test conversion with each coordinate system
            result = spz.ply_to_spz(sample_ply_path, spz_path, coordinate_system=cs)
            
            # With the corrected PLY format, conversion should succeed for all coordinate systems
            assert result is True, f"PLY to SPZ conversion should succeed with {cs} coordinate system"
            
            # Verify the output
            assert os.path.exists(spz_path), f"SPZ file with {cs} coordinate system should be created"
            assert os.path.getsize(spz_path) > 0, f"SPZ file with {cs} coordinate system should not be empty"
            print(f"Conversion with {cs} coordinate system succeeded") 