"""The eclipse fit and the ephemeris fit, on synthetic eclipses and on the windows of the real light curves."""
import numpy as np
import pytest

from hd170312.eclipses import assign_cycles, fit_ephemeris, fit_event_time, predict
from hd170312.events import fit_events


def dip_flux(seed, t0=10.0, depth=0.08, duration=0.75):
    """3 days of 2-minute cadences with white noise and a trapezoid eclipse (ingress half the half-duration) at t0."""
    time = 8.5 + np.arange(0, 3.0, 2 / 1440)
    x = np.abs(time - t0)
    shape = np.clip((duration / 2 - x) / (duration / 4), 0, 1)
    return time, (1 + np.random.default_rng(seed).normal(0, 3e-4, time.size)) * (1 - depth * shape)


# the Note's values: mid-times and errors (BTJD), depth, duration (d)
NOTE = {"S39": (2377.1543, 0.0002, 0.0855, 0.747), "S102": (4131.0663, 0.0001, 0.1536, 0.653), "S103": (4169.8482, 0.0005, 0.0850, 0.741)}


def test_an_injected_eclipse_is_recovered_within_its_error_and_its_depth():
    t0 = 10.03
    time, flux = dip_flux(5, t0=t0)
    found = fit_event_time(time, flux, 10.0)
    assert abs(found["time"] - t0) < 3 * found["error"] and abs(found["depth"] - 0.08) < 0.002 and abs(found["duration"] - 0.75) < 0.02
    assert found["error"] == pytest.approx(np.hypot(found["statistical_error"], found["systematic_error"]))


def test_a_window_with_no_eclipse_is_an_error():
    time, flux = dip_flux(1, depth=0.0)
    with pytest.raises(ValueError, match="no dip"):
        fit_event_time(time, flux, 10.0)


def test_an_eclipse_cut_by_the_end_of_the_data_is_partial_and_is_refused():
    time, flux = dip_flux(2, t0=10.0, depth=0.08, duration=0.75)
    keep = time < 10.05
    with pytest.raises(ValueError, match="partial"):
        fit_event_time(time[keep], flux[keep], 10.0)


def test_the_three_eclipses_of_the_light_curves_are_fitted_as_the_note_gives_them(windows):
    found = fit_events(loader=windows)
    for name, (time, error, depth, duration) in NOTE.items():
        assert found[name]["time"] == pytest.approx(time, abs=5e-5), name
        assert found[name]["error"] == pytest.approx(error, abs=6e-5), name
        assert found[name]["depth"] == pytest.approx(depth, abs=6e-4), name
        assert found[name]["duration"] == pytest.approx(duration, abs=1.5e-3), name
    assert [found[name]["kind"] for name in NOTE] == ["shallow", "deep", "shallow"]


def test_the_ephemeris_recovers_the_period_and_the_epoch_of_each_kind():
    events = [{"time": 100.0, "error": 0.001, "cycle": 0, "kind": "shallow"}, {"time": 300.0, "error": 0.001, "cycle": 4, "kind": "shallow"}, {"time": 125.0, "error": 0.001, "cycle": 0, "kind": "deep"}]
    fit = fit_ephemeris(events)
    assert fit["period"] == pytest.approx(50.0) and fit["epochs"]["deep"] == pytest.approx(125.0) and fit["dof"] == 0 and fit["reduced_chi2"] is None
    assert predict(fit, "deep", 2)[0] == pytest.approx(225.0)


def test_one_event_of_a_kind_cannot_fix_the_period():
    with pytest.raises(ValueError):
        fit_ephemeris([{"time": 1.0, "error": 0.1, "cycle": 0, "kind": "a"}, {"time": 2.0, "error": 0.1, "cycle": 0, "kind": "b"}, {"time": 3.0, "error": 0.1, "cycle": 0, "kind": "b"}])


def test_cycles_are_counted_from_the_first_event_of_each_kind():
    assert assign_cycles([10.0, 12.0, 30.0], ["a", "b", "a"], 10.0) == [0, 0, 2]
