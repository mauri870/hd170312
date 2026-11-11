"""The box search recovers a periodic box dip injected into synthetic sectors."""
import numpy as np
import pytest

pytest.importorskip("transitgpu", reason="uv sync --locked --group search")
from hd170312 import search


def synthetic(sector):
    rng = np.random.default_rng(sector)
    time = 1000.0 + 30.0 * sector + np.arange(0, 20.0, 1 / 72)
    flux = 1 + rng.normal(0, 4e-4, time.size)
    flux[((time - 1001.0) % 3.7) < 0.1] -= 0.004
    return {"time": time, "pdcsap": flux}


def test_the_stack_is_detrended_and_the_search_finds_the_injected_period():
    time, flux = search.stack(loader=synthetic, sectors=(0, 1, 2))
    assert abs(np.median(flux) - 1) < 1e-3 and time.size == flux.size
    found = search.run(time, flux, 2.0, 6.0)
    assert abs(found["period_days"] - 3.7) < 2 * found["grid_step_days"] and found["sde"] > 8 and 0.002e6 < found["depth_ppm"] < 0.006e6

