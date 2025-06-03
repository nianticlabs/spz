# SPZ Python Tests

This directory contains test suites for the SPZ Python bindings, organized by functionality.

## Test Files

### `test_import.py`
Tests for basic module import functionality:
- Module import validation
- Function availability checks
- Version information

### `test_ply_to_spz.py`
Tests for PLY to SPZ conversion:
- File format conversion
- Coordinate system handling
- Error handling with missing files
- Compression verification
- Multiple coordinate systems (RDF, RUB, LUF, RUF)

### `test_spz_to_ply.py`
Tests for SPZ to PLY conversion:
- Reverse conversion functionality
- Coordinate system support
- Error handling
- File validation

### `test_roundtrip.py`
Tests for roundtrip conversions (PLY → SPZ → PLY):
- Data integrity verification
- Compression analysis
- Multiple coordinate systems
- Size ratio validation

### `conftest.py`
Shared test utilities and fixtures:
- `sample_ply_path` fixture for test data
- `create_simple_ply()` utility for creating test PLY files
- Common test helpers

## Running Tests

### All Tests
```bash
# From project root
pytest python/tests/

# With verbose output
pytest python/tests/ -v

# With output capture disabled (to see print statements)
pytest python/tests/ -s
```

### Individual Test Files
```bash
# Test imports only
pytest python/tests/test_import.py

# Test PLY to SPZ conversion
pytest python/tests/test_ply_to_spz.py

# Test SPZ to PLY conversion
pytest python/tests/test_spz_to_ply.py

# Test roundtrip conversions
pytest python/tests/test_roundtrip.py
```

### Specific Tests
```bash
# Run a specific test function
pytest python/tests/test_ply_to_spz.py::test_ply_to_spz_simple_conversion

# Run tests matching a pattern
pytest python/tests/ -k "conversion"
```

## Test Requirements

- Python 3.7+
- pytest
- SPZ Python bindings
- Sample PLY files (in `samples/` directory)

## Test Data

Tests use:
- Sample PLY files from the `samples/` directory
- Dynamically created simple PLY files for basic testing
- Temporary directories for output files (automatically cleaned up)

## Coverage

The tests cover:
- ✅ Basic imports and function availability
- ✅ PLY to SPZ conversion with various coordinate systems
- ✅ SPZ to PLY conversion with various coordinate systems
- ✅ Roundtrip conversion integrity
- ✅ Error handling for invalid inputs
- ✅ File compression verification
- ✅ Coordinate system validation 