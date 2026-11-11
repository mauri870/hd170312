"""The lower bound on the radius of the body B that eclipses A: with and without limb darkening, for the radii of A of the TIC and of the two-star fits."""
import pytest

from hd170312 import radius
from hd170312.events import fit_events


def test_the_limb_darkening_is_interpolated_in_teff_between_the_two_rows_of_claret_2017_and_not_extrapolated():
    a, b = radius.limb_darkening(6645.0)
    fraction = (6645.0 - 6500.0) / 250.0
    assert a == pytest.approx(0.2028 + fraction * (0.1843 - 0.2028)) and b == pytest.approx(0.3226 + fraction * (0.3280 - 0.3226))
    assert radius.limb_darkening(6500.0) == pytest.approx((0.2028, 0.3226)) and radius.limb_darkening(6750.0) == pytest.approx((0.1843, 0.3280))
    with pytest.raises(ValueError, match="6500"):
        radius.limb_darkening(6400.0)


def test_the_centre_is_brighter_than_the_mean_of_the_disk_by_one_over_one_minus_a_third_of_a_and_a_sixth_of_b():
    assert radius.central_over_mean(0.0, 0.0) == 1.0
    assert radius.central_over_mean(0.6, 0.0) == pytest.approx(1 / 0.8)
    assert radius.central_over_mean(0.192, 0.3257) == pytest.approx(1.1342, abs=1e-3)


def test_the_bound_is_the_square_root_of_the_depth_over_the_brightness_of_the_centre_times_the_radius_of_a():
    assert radius.lower_bound(0.04, 2.0) == pytest.approx(0.4)
    assert radius.lower_bound(0.04, 2.0, (0.6, 0.0)) == pytest.approx(0.4 * 0.8**0.5), "limb darkening lowers the bound: the centre is brighter than the mean"
    with pytest.raises(ValueError):
        radius.lower_bound(-0.01, 2.0)


def test_the_table_has_each_kind_each_radius_with_and_without_limb_darkening_and_the_conservative_case_is_the_smallest(windows):
    table = radius.table(fit_events(loader=windows))
    assert len(table["cases"]) == 12 and {case["kind"] for case in table["cases"]} == {"shallow", "deep"}
    assert {case["radius_of_a"] for case in table["cases"]} == {3.46, 2.64, 2.27}
    for case in table["cases"]:
        assert case["radius_of_b"] > 0 and case["limb_darkening"] in (True, False)
    smallest = min(table["cases"], key=lambda case: case["radius_of_b"])
    assert smallest["kind"] == "shallow" and smallest["radius_of_a"] == 2.27 and smallest["limb_darkening"] is True
    assert table["conservative"]["radius_of_b"] == smallest["radius_of_b"] == pytest.approx(0.62, abs=0.01)
    shallow = next(case for case in table["cases"] if case["kind"] == "shallow")["depth"]
    assert shallow == pytest.approx(0.0850, abs=1e-3), "the shallower of the two shallow eclipses"
