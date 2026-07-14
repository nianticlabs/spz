"""Tests for PackOptions.fractional_bits."""

import struct

import numpy as np
import pytest

import spz


def _cloud(num_points=1000, extent=8.0, seed=0):
    rng = np.random.default_rng(seed)
    c = spz.GaussianCloud()
    c.sh_degree = 1
    c.positions = rng.uniform(-extent, extent, size=num_points * 3).astype(np.float32)
    c.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    c.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    c.alphas = rng.uniform(-3.0, 3.0, size=num_points).astype(np.float32)
    c.colors = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    c.sh = rng.uniform(-0.5, 0.5, size=num_points * 9).astype(np.float32)
    return c


def _save(cloud, tmp_path, **opts):
    o = spz.PackOptions()
    for k, v in opts.items():
        setattr(o, k, v)
    _save.counter += 1
    path = str(tmp_path / f"out_{_save.counter}.spz")
    return spz.save_spz(cloud, o, path), path

_save.counter = 0


def test_defaults_and_setters():
    assert spz.MIN_FRACTIONAL_BITS == 4
    assert spz.MAX_FRACTIONAL_BITS == 23
    assert spz.DEFAULT_FRACTIONAL_BITS == 12
    o = spz.PackOptions()
    assert o.fractional_bits == 12
    o.fractional_bits = 8
    assert o.fractional_bits == 8


@pytest.mark.parametrize("n", [4, 8, 12, 16, 20, 23])
def test_fractional_bits_round_trip(n, tmp_path):
    """Header echoes fractional_bits, and position error <= 2^-N."""
    # Pick extent that fits in +/- 2^(23-N) (the representable range), capped
    # at 8.0 for ergonomics. Without this cap, the encoder would refuse to
    # save at N=23 (range +/- 1.0).
    extent = min(8.0, 0.5 * 2.0 ** (23 - n))
    cloud = _cloud(extent=extent)
    orig = np.asarray(cloud.positions, dtype=np.float32).copy()
    ok, path = _save(cloud, tmp_path, fractional_bits=n)
    assert ok
    # NgspFileHeader layout: magic[0..3], version[4..7], numPoints[8..11],
    # shDegree[12], fractionalBits[13]. fractionalBits should echo `n`.
    with open(path, "rb") as f:
        magic, version, _, _, frac_bits = struct.unpack("<IIIBB", f.read(14))
    assert magic == 0x5053474E and version == 4 and frac_bits == n
    # Round-trip error bound: rounding to nearest at step 2^-N.
    recovered = np.asarray(spz.load_spz(path).positions, dtype=np.float32)
    max_err = float(np.max(np.abs(recovered - orig)))
    assert max_err <= 2.0 ** -n, f"max_err={max_err:.6g} > 2^-{n}"


@pytest.mark.parametrize("n", [
    spz.MIN_FRACTIONAL_BITS - 1,
    0,
    spz.MAX_FRACTIONAL_BITS + 1,
])
def test_fractional_bits_out_of_range_rejected(n, tmp_path):
    cloud = _cloud(num_points=10)
    ok, _ = _save(cloud, tmp_path, fractional_bits=n)
    assert ok is False
