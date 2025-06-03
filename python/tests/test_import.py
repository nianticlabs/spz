"""Tests for basic SPZ module import functionality."""

import pytest
import spz


def test_import():
    """Test that the module can be imported and has expected attributes."""
    assert hasattr(spz, 'ply_to_spz')
    assert hasattr(spz, 'spz_to_ply')
    assert hasattr(spz, '__version__')


def test_import_ply_to_spz():
    """Test that the ply_to_spz function is available."""
    assert hasattr(spz, 'ply_to_spz')
    assert callable(spz.ply_to_spz)


def test_import_spz_to_ply():
    """Test that the spz_to_ply function is available."""
    assert hasattr(spz, 'spz_to_ply')
    assert callable(spz.spz_to_ply) 