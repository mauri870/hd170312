"""The alias rule on a synthetic star whose eclipses repeat every 30 d."""
import numpy as np
import pytest

from hd170312 import alias

CADENCE = 2 / 1440


def synthetic_star():
    """Three sectors of 27 d with 300 ppm noise and an eclipse every 30 d: a deep one (15 percent) at 15 + 30 n and a shallow one (8 percent) at 22.5 + 30 n."""
    generator = np.random.default_rng(3)
    sectors = {}
    for number, start in enumerate((0.0, 30.0, 60.0)):
        time = np.arange(start, start + 27.0, CADENCE)
        flux = 1 + generator.normal(0, 300e-6, time.size)
        for center, depth in ((start + 15.0, 0.15), (start + 22.5, 0.08)):
            x = np.abs(time - center)
            flux = flux * (1 - depth * np.clip((0.35 - x) / 0.2, 0, 1))
        sectors[number + 1] = (time, flux)
    return sectors


def test_the_period_of_the_eclipses_and_its_half_survive_and_the_others_are_rejected_by_an_eclipse_that_would_have_been_seen():
    sectors = synthetic_star()
    deepest = {"deep": [(15.0, 1)], "shallow": [(22.5, 1), (52.5, 2)]}
    prepared = alias.prepare(sectors, deepest, [15.0, 22.5, 45.0, 52.5, 75.0, 82.5])
    table = alias.alias_table(sectors, prepared, anchors=(("deep", 15.0), ("shallow", 82.5)), separation=60.0, n_values=range(1, 7), event_times=[15.0, 22.5, 45.0, 52.5, 75.0, 82.5])
    assert alias.survivors(table) == [1, 2]
    assert table[3]["rejected"] and table[3]["refuting_eclipses"][0]["share_of_box"] >= 0.5
    assert table[2]["tested"] >= 4 and table[2]["of_which_events"] == table[2]["eclipses_in_the_data"]


def test_an_eclipse_with_less_than_half_its_box_covered_is_not_tested():
    time = np.arange(0.0, 5.0, CADENCE)
    flux = np.ones_like(time)
    prepared = {"grid": np.array([-1.0, 1.0]), "profile": {"deep": np.array([0.1, 0.1])}, "observed_depth": {"deep": 0.1}, "sigma_box": {1: 1e-4}}
    record = alias.test_eclipse({1: (time, flux)}, prepared, "deep", 5.1, 1)
    assert not record["tested"] and record["share_of_box"] < 0.5


def test_a_missing_dip_that_an_injection_would_have_shown_is_a_rejection_and_a_noisy_one_is_not():
    time = np.arange(0.0, 10.0, CADENCE)
    generator = np.random.default_rng(1)
    prepared = {"grid": np.array([-30.0, 30.0]), "profile": {"deep": np.array([0.1, 0.1])}, "observed_depth": {"deep": 0.1}, "sigma_box": {1: 1e-3}}
    flat = alias.test_eclipse({1: (time, 1 + generator.normal(0, 3e-4, time.size))}, prepared, "deep", 5.0, 1)
    assert flat["tested"] and not flat["present"] and flat["recovery_detected"]
    dip = 1 + generator.normal(0, 3e-4, time.size) - 0.1 * (np.abs(time - 5.0) < 0.3)
    present = alias.test_eclipse({1: (time, dip)}, prepared, "deep", 5.0, 1)
    assert present["present"]


def test_a_period_that_does_not_predict_the_fourth_eclipse_is_rejected_by_it_and_the_two_that_do_are_not():
    anchors = (("deep", alias.T_DEEP), ("shallow", alias.T_SHALLOW_2))
    separation = alias.T_SHALLOW_2 - alias.T_SHALLOW_1
    table = alias.fourth_event_table(range(1, 31), anchors, separation, "deep", 3855.2675, 0.001)
    assert [n for n, row in table.items() if row["fits"]] == [13, 26]
    assert table[26]["residual_days"] == pytest.approx(0.0, abs=5e-4) and table[26]["cycles"] == -4 and table[13]["cycles"] == -2
    assert abs(table[1]["residual_days"]) == pytest.approx(275.8, abs=0.1) and not table[9]["fits"]
    assert all(row["fits"] for row in alias.fourth_event_table([1, 2, 9], anchors, separation, "deep", 3855.2675, 400.0).values()), "an error of the event of the size of the residual accepts every period"


def test_the_observed_depth_is_the_mean_of_the_profile_within_half_the_box_in_hours():
    sectors = synthetic_star()
    deepest = {"deep": [(15.0, 1)], "shallow": [(22.5, 1), (52.5, 2)]}
    prepared = alias.prepare(sectors, deepest, [15.0, 22.5, 45.0, 52.5, 75.0, 82.5])
    assert alias.BOX_HOURS == 13.45 and alias.BOX == pytest.approx(alias.BOX_HOURS / 24)
    for kind, profile in prepared["profile"].items():
        assert prepared["observed_depth"][kind] == pytest.approx(profile[np.abs(prepared["grid"]) < alias.BOX_HOURS / 2].mean())
