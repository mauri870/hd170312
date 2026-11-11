"""The eclipse constants written into the analysis modules are the fit's output to the digits written."""
import numpy as np
import pytest

from hd170312 import alias, figure, pooled, preregistered
from hd170312.events import fit_events


@pytest.fixture(scope="module")
def fit(windows):
    return fit_events(loader=windows)


def test_the_mid_times_in_the_modules_are_the_fitted_times_to_the_digits_written(fit):
    s39, s102, s103 = fit["S39"]["time"], fit["S102"]["time"], fit["S103"]["time"]
    assert alias.T_SHALLOW_1 == pytest.approx(s39, abs=5e-5) and alias.T_DEEP == pytest.approx(s102, abs=5e-5) and alias.T_SHALLOW_2 == pytest.approx(s103, abs=5e-5)
    assert (preregistered.KINDS[0][0], preregistered.KINDS[1][0]) == (alias.T_SHALLOW_1, alias.T_DEEP)
    assert (figure.DEEP_ANCHOR, figure.SHALLOW_ANCHOR) == (alias.T_DEEP, alias.T_SHALLOW_2)
    assert figure.DELTA_T == pytest.approx(s103 - s39, abs=5e-4) and preregistered.P26 * 26 == pytest.approx(figure.DELTA_T, abs=1e-9)


def test_the_depths_durations_and_ingress_in_the_modules_are_the_fitted_ones_to_the_digits_written(fit):
    shallow, deep = fit["S39"], fit["S102"]
    assert preregistered.KINDS[0][1] == pytest.approx(shallow["depth"], abs=5e-5) and preregistered.KINDS[1][1] == pytest.approx(deep["depth"], abs=5e-5)
    assert preregistered.KINDS[0][2] == pytest.approx(shallow["duration"], abs=5e-4) and preregistered.KINDS[1][2] == pytest.approx(deep["duration"], abs=5e-4)
    assert pooled.SHAPES["shallow"][0] == preregistered.KINDS[0][1] and pooled.SHAPES["deep"][0] == preregistered.KINDS[1][1]
    assert pooled.SHAPES["shallow"][1] == preregistered.KINDS[0][2] and pooled.SHAPES["deep"][1] == preregistered.KINDS[1][2]
    assert pooled.SHAPES["shallow"][2] == pytest.approx(shallow["ingress_fraction"], abs=5e-4) and pooled.SHAPES["deep"][2] == pytest.approx(deep["ingress_fraction"], abs=5e-4)


def test_the_errors_of_the_anchors_in_the_alias_test_are_the_fitted_errors_to_the_digits_written(fit):
    assert alias.SIGMA_ANCHOR["deep"] == pytest.approx(fit["S102"]["error"], abs=5e-5)
    assert alias.SIGMA_ANCHOR["shallow"] == pytest.approx(fit["S103"]["error"], abs=5e-5)
    assert alias.SIGMA_SEPARATION == pytest.approx(np.hypot(fit["S39"]["error"], fit["S103"]["error"]), abs=5e-5)


def test_the_deepest_points_of_the_events_in_the_alias_test_are_within_a_few_minutes_of_the_fitted_times(fit):
    deepest = {102: alias.DEEPEST["deep"][0][0], 39: alias.DEEPEST["shallow"][0][0], 103: alias.DEEPEST["shallow"][1][0]}
    fitted = {102: fit["S102"]["time"], 39: fit["S39"]["time"], 103: fit["S103"]["time"]}
    for sector, time in deepest.items():
        assert abs(time - fitted[sector]) < 0.02, (sector, time, fitted[sector])
