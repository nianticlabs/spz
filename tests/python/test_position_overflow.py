"""save_spz must reject positions that overflow the int24 fixed-point range.

Positions are stored as 24-bit signed fixed-point with fractionalBits=12 by
default, so the per-axis representable range is +/- 2^(23 - 12) = 2048 metres.
Values outside that range silently wrap (low 24 bits, sign-extended) on
existing code. This test asserts the encoder now refuses the save instead.
"""
import os

import numpy as np
import pytest

import spz

# The encoder uses fractionalBits = 12. int24 range is [-2^23, 2^23 - 1]
# so the largest representable input is (2^23 - 1) / 2^12 = 2047.999755859375.
# Anything strictly larger overflows after rounding.
MAX_FIXED = (1 << 23) - 1
MIN_FIXED = -(1 << 23)
SCALE = 1 << 12
MAX_SAFE = MAX_FIXED / SCALE       # 2047.999755859375
MIN_SAFE = MIN_FIXED / SCALE       # -2048.0


def _cloud_with_position(x, y=0.0, z=0.0):
    """Single-point cloud with the given position; everything else zero."""
    c = spz.GaussianCloud()
    c.positions = np.array([x, y, z], dtype=np.float32)
    c.scales = np.zeros(3, dtype=np.float32)
    c.rotations = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float32)
    c.alphas = np.zeros(1, dtype=np.float32)
    c.colors = np.zeros(3, dtype=np.float32)
    return c


@pytest.mark.parametrize("v", [
    MAX_SAFE,                  # exactly the largest representable value
    MAX_SAFE - 1.0,            # comfortably inside
    0.0,                       # zero
    -2048.0,                   # min int24 / 2^12 (signed range is asymmetric)
    -100.0,
])
def test_in_range_positions_save_successfully(v, tmp_path):
    cloud = _cloud_with_position(v)
    path = str(tmp_path / "ok.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), path) is True
    # Round-trip preserves the value within the quantization step.
    loaded = spz.load_spz(path)
    recovered = float(np.asarray(loaded.positions)[0])
    assert abs(recovered - v) <= 1.0 / SCALE


@pytest.mark.parametrize("v", [
    2048.0,                    # smallest float that overflows int24 after rounding
    2050.0,                    # comfortably overflowing
    -2048.001,                 # just past the negative limit (rounds to -8388612)
    -10000.0,                  # heavy overflow on the negative side
    1e6,                       # very large
    -1e6,
])
def test_out_of_range_positions_are_rejected(v, tmp_path):
    cloud = _cloud_with_position(v)
    path = str(tmp_path / "fail.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), path) is False
    # No file on disk on failure.
    assert not os.path.exists(path) or os.path.getsize(path) == 0


def test_rounding_edge_case_rejected(tmp_path):
    """A value just below the boundary that *rounds up* into overflow must be rejected.

    With fractionalBits=12, scale=4096, the integer ceiling is 8388607. Any
    float V where round(V * 4096) >= 8388608 overflows; the smallest such V is
    8388607.5 / 4096 = 2047.9998779296875. The test value below sits right at
    that boundary -- a naive |V| < 2048 check would miss it.
    """
    v = 8388607.5 / SCALE        # rounds to 8388608 (overflows by 1)
    cloud = _cloud_with_position(v)
    path = str(tmp_path / "edge.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), path) is False


def test_existing_samples_still_save(tmp_path):
    """The two repo samples are well within range; they must continue to save."""
    for name in ["hornedlizard.spz", "racoonfamily.spz"]:
        sample = f"samples/{name}"
        if not os.path.exists(sample):
            pytest.skip(f"sample missing: {sample}")
        cloud = spz.load_spz(sample)
        out = str(tmp_path / name)
        assert spz.save_spz(cloud, spz.PackOptions(), out) is True
        assert os.path.getsize(out) > 0
