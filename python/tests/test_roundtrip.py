"""Tests for roundtrip conversion functionality (PLY -> SPZ -> PLY)."""

import pytest
import os
import tempfile
import spz
from .conftest import create_simple_ply


def test_spz_to_ply_roundtrip(sample_ply_path):
    """Test roundtrip conversion: PLY -> SPZ -> PLY."""
    # Skip test if sample file doesn't exist
    if not os.path.exists(sample_ply_path):
        pytest.skip("Sample PLY file not found")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        spz_path = os.path.join(tmpdir, "intermediate.spz")
        output_ply_path = os.path.join(tmpdir, "roundtrip.ply")
        
        # PLY -> SPZ -> PLY roundtrip
        result1 = spz.ply_to_spz(sample_ply_path, spz_path, coordinate_system="RDF")
        assert result1, "PLY to SPZ conversion should succeed"
        
        result2 = spz.spz_to_ply(spz_path, output_ply_path, coordinate_system="RDF")
        assert result2, "SPZ to PLY conversion should succeed"
        
        # Both files should exist and be non-empty
        assert os.path.exists(spz_path), "SPZ file should exist"
        assert os.path.exists(output_ply_path), "Output PLY file should exist"
        assert os.path.getsize(spz_path) > 0, "SPZ file should not be empty"
        assert os.path.getsize(output_ply_path) > 0, "Output PLY file should not be empty"
        
        # Check file sizes
        original_size = os.path.getsize(sample_ply_path)
        spz_size = os.path.getsize(spz_path)
        output_size = os.path.getsize(output_ply_path)
        
        print(f"Roundtrip test:")
        print(f"  Original PLY: {original_size:,} bytes")
        print(f"  SPZ: {spz_size:,} bytes")
        print(f"  Roundtrip PLY: {output_size:,} bytes")
        print(f"  Compression ratio: {original_size/spz_size:.2f}x")
        
        # Output PLY might be slightly different size due to precision changes,
        # but should be in the same ballpark
        size_ratio = output_size / original_size
        assert 0.8 <= size_ratio <= 1.2, f"Roundtrip PLY size should be similar to original (ratio: {size_ratio:.2f})"


def test_simple_roundtrip():
    """Test roundtrip with a simple created PLY file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        original_ply_path = os.path.join(tmpdir, "original.ply")
        spz_path = os.path.join(tmpdir, "intermediate.spz")
        roundtrip_ply_path = os.path.join(tmpdir, "roundtrip.ply")
        
        # Create original PLY file
        create_simple_ply(original_ply_path)
        
        # Roundtrip conversion
        result1 = spz.ply_to_spz(original_ply_path, spz_path, coordinate_system="RDF")
        assert result1, "PLY to SPZ conversion should succeed"
        
        result2 = spz.spz_to_ply(spz_path, roundtrip_ply_path, coordinate_system="RDF")
        assert result2, "SPZ to PLY conversion should succeed"
        
        # Verify all files exist and have content
        assert os.path.exists(spz_path), "SPZ file should exist"
        assert os.path.exists(roundtrip_ply_path), "Roundtrip PLY file should exist"
        assert os.path.getsize(spz_path) > 0, "SPZ file should not be empty"
        assert os.path.getsize(roundtrip_ply_path) > 0, "Roundtrip PLY file should not be empty"
        
        # Check compression
        original_size = os.path.getsize(original_ply_path)
        spz_size = os.path.getsize(spz_path)
        roundtrip_size = os.path.getsize(roundtrip_ply_path)
        
        print(f"Simple roundtrip test:")
        print(f"  Original PLY: {original_size:,} bytes")
        print(f"  SPZ: {spz_size:,} bytes")
        print(f"  Roundtrip PLY: {roundtrip_size:,} bytes")
        print(f"  Compression ratio: {original_size/spz_size:.2f}x")
        
        assert spz_size < original_size, "SPZ should be smaller than original PLY"


def test_roundtrip_different_coordinate_systems(sample_ply_path):
    """Test roundtrip conversions with different coordinate systems."""
    if not os.path.exists(sample_ply_path):
        pytest.skip("Sample PLY file not found")
    
    coordinate_systems = ["RDF", "RUB", "LUF", "RUF"]
    
    with tempfile.TemporaryDirectory() as tmpdir:
        for cs in coordinate_systems:
            spz_path = os.path.join(tmpdir, f"intermediate_{cs}.spz")
            output_ply_path = os.path.join(tmpdir, f"roundtrip_{cs}.ply")
            
            # Roundtrip with specific coordinate system
            result1 = spz.ply_to_spz(sample_ply_path, spz_path, coordinate_system=cs)
            assert result1, f"PLY to SPZ conversion with {cs} should succeed"
            
            result2 = spz.spz_to_ply(spz_path, output_ply_path, coordinate_system=cs)
            assert result2, f"SPZ to PLY conversion with {cs} should succeed"
            
            # Verify files
            assert os.path.exists(spz_path), f"SPZ file with {cs} should exist"
            assert os.path.exists(output_ply_path), f"Output PLY file with {cs} should exist"
            assert os.path.getsize(spz_path) > 0, f"SPZ file with {cs} should not be empty"
            assert os.path.getsize(output_ply_path) > 0, f"Output PLY file with {cs} should not be empty"
            
            print(f"Roundtrip with {cs} coordinate system succeeded") 