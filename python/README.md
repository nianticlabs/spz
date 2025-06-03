# SPZ Python Package

This directory contains the Python components of the SPZ (compressed Gaussian splat) format library.

## Directory Structure

```
python/
├── spz/                    # Main Python package
│   ├── __init__.py        # Package initialization and public API
│   ├── core.py           # Core Python utilities
│   ├── _core.cpython-*.so # Compiled C++ extension
│   └── bindings/         # C++ binding source code
│       └── spz_bindings.cpp
├── examples/             # Example scripts
│   ├── README.md        # Examples documentation
│   ├── example.py       # Basic PLY to SPZ conversion
│   └── example_bidirectional.py # Advanced bidirectional conversion
└── tests/               # Test suites
    ├── README.md       # Test documentation
    ├── conftest.py     # Shared test utilities
    ├── test_import.py  # Import functionality tests
    ├── test_ply_to_spz.py   # PLY to SPZ conversion tests
    ├── test_spz_to_ply.py   # SPZ to PLY conversion tests
    └── test_roundtrip.py    # Roundtrip conversion tests
```

## Package Components

### Main Package (`spz/`)
- **Core module**: Provides the main Python API for SPZ format conversion
- **Bindings**: C++ extension providing high-performance conversion routines
- **Public API**: 
  - `spz.ply_to_spz()` - Convert PLY files to compressed SPZ format
  - `spz.spz_to_ply()` - Convert SPZ files back to PLY format

### Examples (`examples/`)
Ready-to-run example scripts demonstrating:
- Basic file conversion
- Bidirectional conversion workflows
- Performance measurement
- Error handling

### Tests (`tests/`)
Comprehensive test suite covering:
- Module imports and API availability
- Conversion functionality for both directions
- Roundtrip data integrity
- Error handling and edge cases
- Multiple coordinate system support

## Installation

From the project root:

```bash
# Install in development mode
pip install -e .

# Or install from wheel
pip install .
```

## Quick Start

```python
import spz

# Convert PLY to SPZ
success = spz.ply_to_spz("input.ply", "output.spz", coordinate_system="RDF")

# Convert SPZ back to PLY
success = spz.spz_to_ply("input.spz", "output.ply", coordinate_system="RDF")
```

## Coordinate Systems

Supported coordinate systems:
- **RDF**: Right-Down-Forward (typical for PLY files)
- **RUB**: Right-Up-Back
- **LUF**: Left-Up-Forward  
- **RUF**: Right-Up-Forward
- **UNSPECIFIED**: Use default handling

## Development

### Building
The C++ extension is built automatically during installation using pybind11.

### Testing
```bash
# Run all tests
pytest python/tests/

# Run with coverage
pytest python/tests/ --cov=spz

# Run specific test categories
pytest python/tests/test_ply_to_spz.py
```

### Examples
```bash
# Basic conversion
python python/examples/example.py samples/splat.ply output.spz

# Advanced features
python python/examples/example_bidirectional.py samples/splat.ply output/
``` 