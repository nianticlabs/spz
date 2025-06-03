# SPZ Python Examples

This directory contains example scripts demonstrating how to use the SPZ Python bindings.

## Examples

### `example.py`
Basic example showing how to convert a PLY file to SPZ format.

**Usage:**
```bash
python example.py <input.ply> <output.spz>
```

**Features:**
- Converts PLY to SPZ using RDF coordinate system
- Shows compression statistics
- Error handling for missing files

### `example_bidirectional.py`
Advanced example demonstrating bidirectional conversion between PLY and SPZ formats.

**Usage:**
```bash
python example_bidirectional.py <input.ply> <output_directory>
```

**Features:**
- PLY → SPZ → PLY roundtrip conversion
- Multiple coordinate system support
- Detailed compression analysis
- Performance timing
- Quality verification

## Running Examples

Make sure you have the SPZ package installed:

```bash
pip install -e .
```

Then you can run the examples from the project root:

```bash
# Basic conversion
python python/examples/example.py samples/splat.ply output.spz

# Bidirectional conversion
python python/examples/example_bidirectional.py samples/splat.ply output/
```

## Requirements

- Python 3.7+
- SPZ Python bindings
- Sample PLY files (available in the `samples/` directory) 