"""SPZ (Sparse Point-cloud Zip) Python bindings.

This module provides Python bindings for the SPZ C++ library,
which handles sparse point-cloud compression and decompression.
"""

try:
    from .spz_bindings import *  # type: ignore[misc]
except ImportError as e:
    # The C++ extension might not be built yet
    import warnings
    warnings.warn(f"SPZ bindings not available: {e}")
