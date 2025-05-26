#!/usr/bin/env python3
"""
Example usage of the SPZ Python bindings.

This example shows how to convert a PLY file to SPZ format.
"""

import sys
import os
import spz

def main():
    if len(sys.argv) != 3:
        print("Usage: python example.py <input.ply> <output.spz>")
        sys.exit(1)
    
    ply_path = sys.argv[1]
    spz_path = sys.argv[2]
    
    if not os.path.exists(ply_path):
        print(f"Error: Input file {ply_path} does not exist")
        sys.exit(1)
    
    print(f"Converting {ply_path} to {spz_path}...")
    
    # Convert PLY to SPZ using RDF coordinate system (typical for PLY files)
    success = spz.ply_to_spz(ply_path, spz_path, coordinate_system="RDF")
    
    if success:
        print("Conversion successful!")
        
        # Check file sizes
        ply_size = os.path.getsize(ply_path)
        spz_size = os.path.getsize(spz_path)
        compression_ratio = ply_size / spz_size if spz_size > 0 else float('inf')
        
        print(f"PLY file size: {ply_size:,} bytes")
        print(f"SPZ file size: {spz_size:,} bytes")
        print(f"Compression ratio: {compression_ratio:.2f}x")
    else:
        print("Conversion failed!")
        sys.exit(1)

if __name__ == "__main__":
    main() 