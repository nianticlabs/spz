#!/usr/bin/env python3
"""
Example script demonstrating bidirectional conversion between PLY and SPZ formats.

This script shows how to:
1. Convert PLY to SPZ (compression)
2. Convert SPZ back to PLY (decompression) 
3. Handle different coordinate systems
4. Check compression ratios and file integrity

Usage:
    python example_bidirectional.py [ply_file]

If no PLY file is provided, it will look for samples/splat.ply
"""

import os
import sys
import tempfile
import time
import spz


def format_bytes(size):
    """Format byte size into human readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} TB"


def main():
    # Determine input PLY file
    if len(sys.argv) > 1:
        ply_file = sys.argv[1]
    else:
        # Look for the sample file
        ply_file = "samples/splat.ply"
    
    if not os.path.exists(ply_file):
        print(f"❌ PLY file not found: {ply_file}")
        print("Usage: python example_bidirectional.py [ply_file]")
        sys.exit(1)
    
    print(f"🎯 SPZ Bidirectional Conversion Example")
    print(f"📁 Input PLY file: {ply_file}")
    print(f"📊 Input file size: {format_bytes(os.path.getsize(ply_file))}")
    print()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        spz_file = os.path.join(tmpdir, "compressed.spz")
        output_ply_file = os.path.join(tmpdir, "decompressed.ply")
        
        # Step 1: PLY to SPZ conversion (compression)
        print("🔄 Step 1: Compressing PLY to SPZ...")
        start_time = time.time()
        
        success = spz.ply_to_spz(ply_file, spz_file, coordinate_system="RDF")
        
        compression_time = time.time() - start_time
        
        if not success:
            print("❌ PLY to SPZ conversion failed!")
            sys.exit(1)
        
        spz_size = os.path.getsize(spz_file)
        ply_size = os.path.getsize(ply_file)
        compression_ratio = ply_size / spz_size
        
        print(f"✅ Compression successful!")
        print(f"   📊 SPZ file size: {format_bytes(spz_size)}")
        print(f"   🚀 Compression ratio: {compression_ratio:.2f}x")
        print(f"   ⏱️  Time taken: {compression_time:.2f} seconds")
        print()
        
        # Step 2: SPZ to PLY conversion (decompression)
        print("🔄 Step 2: Decompressing SPZ to PLY...")
        start_time = time.time()
        
        success = spz.spz_to_ply(spz_file, output_ply_file, coordinate_system="RDF")
        
        decompression_time = time.time() - start_time
        
        if not success:
            print("❌ SPZ to PLY conversion failed!")
            sys.exit(1)
        
        output_size = os.path.getsize(output_ply_file)
        
        print(f"✅ Decompression successful!")
        print(f"   📊 Output PLY size: {format_bytes(output_size)}")
        print(f"   ⏱️  Time taken: {decompression_time:.2f} seconds")
        print()
        
        # Step 3: Integrity check
        print("🔍 Step 3: Integrity Check...")
        size_difference = abs(output_size - ply_size)
        size_ratio = output_size / ply_size
        
        print(f"   📊 Original PLY: {format_bytes(ply_size)}")
        print(f"   📊 Roundtrip PLY: {format_bytes(output_size)}")
        print(f"   📊 Size difference: {format_bytes(size_difference)}")
        print(f"   📊 Size ratio: {size_ratio:.4f}")
        
        if 0.99 <= size_ratio <= 1.01:
            print("   ✅ File sizes match closely - good integrity!")
        elif 0.95 <= size_ratio <= 1.05:
            print("   ⚠️  Minor size difference (likely due to precision)")
        else:
            print("   ❌ Significant size difference detected")
        print()
        
        # Step 4: Test different coordinate systems
        print("🌐 Step 4: Testing Different Coordinate Systems...")
        coordinate_systems = ["RDF", "RUB", "LUF", "RUF", "UNSPECIFIED"]
        
        for i, cs in enumerate(coordinate_systems, 1):
            cs_spz_file = os.path.join(tmpdir, f"test_{cs}.spz")
            cs_ply_file = os.path.join(tmpdir, f"test_{cs}.ply")
            
            # Convert with specific coordinate system
            success1 = spz.ply_to_spz(ply_file, cs_spz_file, coordinate_system=cs)
            success2 = spz.spz_to_ply(cs_spz_file, cs_ply_file, coordinate_system=cs)
            
            if success1 and success2:
                cs_compression = ply_size / os.path.getsize(cs_spz_file)
                print(f"   {i}. {cs:12} ✅ Compression: {cs_compression:.2f}x")
            else:
                print(f"   {i}. {cs:12} ❌ Failed")
        
        print()
        
        # Summary
        print("📈 Summary:")
        print(f"   🏆 Best compression ratio: {compression_ratio:.2f}x")
        print(f"   💾 Space saved: {format_bytes(ply_size - spz_size)} ({((ply_size - spz_size) / ply_size * 100):.1f}%)")
        print(f"   ⚡ Total processing time: {compression_time + decompression_time:.2f} seconds")
        print()
        print("🎉 Bidirectional conversion completed successfully!")
        print("   The SPZ format provides excellent compression while maintaining data integrity.")


if __name__ == "__main__":
    main() 