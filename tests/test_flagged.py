"""The eclipse of sector 93 and the egress of sector 101 in the SAP flux, on synthetic eclipses and on the windows of the real light curves."""
from pathlib import Path

import numpy as np
import pytest

from hd170312 import flagged
from hd170312.events import fit_events

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="module")
def sap_loader(windows):
    """sector -> {time, sap, quality}: the complete eclipses from the fixture windows (QUALITY 0), sectors 93 and 101 from the fixture with the flagged cadences."""
    arrays = np.load(FIXTURES / "flagged.npz")

    def load(sector):
        if sector in (93, 101):
            return {"time": arrays[f"time_{sector}"], "sap": arrays[f"sap_{sector}"], "quality": arrays[f"quality_{sector}"]}
        columns = windows(sector)
        return {"time": columns["time"], "sap": columns["sap"], "quality": np.zeros(columns["time"].size, dtype=int)}

    return load


def trapezoid_flux(t0, depth, duration=0.75, fraction=0.8, start=8.5, length=3.0, seed=0):
    time = start + np.arange(0, length, 2 / 1440)
    shape = np.clip((duration / 2 - np.abs(time - t0)) / (fraction * duration / 2), 0, 1)
    return time, (1 + np.random.default_rng(seed).normal(0, 3e-4, time.size)) * (1 - depth * shape)


def test_a_template_with_a_fixed_shape_recovers_the_time_and_the_depth_of_an_injected_eclipse():
    time, flux = trapezoid_flux(10.02, 0.08)
    found = flagged.fit_template(time, flux, (0.75, 0.8), centre=10.0, start=(10.0, 0.05))
    assert found["time"] == pytest.approx(10.02, abs=2e-3) and found["depth"] == pytest.approx(0.08, abs=1e-3)
    fixed = flagged.fit_template(time, flux, (0.75, 0.8), centre=10.0, start=(10.02, 0.05), fixed_time=10.02)
    assert fixed["time"] == 10.02 and fixed["depth"] == pytest.approx(0.08, abs=1e-3)


def test_a_window_with_fewer_than_a_hundred_cadences_is_an_error():
    time, flux = trapezoid_flux(10.0, 0.08, length=0.1)
    with pytest.raises(ValueError, match="cadences"):
        flagged.fit_template(time, flux, (0.75, 0.8), centre=10.0, start=(10.0, 0.05))


def test_the_flagged_cadences_of_sector_93_hold_the_deep_eclipse_at_the_predicted_time_and_depth(sap_loader):
    found = flagged.fourth_eclipse(sap_loader)
    assert found["time"] == pytest.approx(3855.2675, abs=5e-4) and found["depth_ratio"] == pytest.approx(1.003, abs=5e-3)
    assert found["flagged_cadences"] > 300, "a third of the fitted cadences are flagged for scattered light"
    assert found["template"]["depth"] == pytest.approx(0.1548, abs=1e-3)


def test_the_egress_of_sector_101_has_the_shallow_depth_at_the_time_that_only_68_95_d_predicts(sap_loader, windows):
    events = fit_events(loader=windows)
    found = flagged.egress(sap_loader, events)
    assert found["predicted_time"] == pytest.approx(4100.8984, abs=1e-4)
    assert found["depth_ratio"] == pytest.approx(0.997, abs=5e-3) and found["gain_over_flat"] > 0.98
    assert found["first_cadence"] == pytest.approx(4101.112, abs=1e-3), "the sector's first finite SAP cadence: only the egress is in the data"
