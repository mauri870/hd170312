"""The NEOWISE-R L1b reading of HD 170312 on synthetic IRSA results: the cuts, the visits, the time, the numbers of the visit with the window and the comparison stars."""
import numpy as np
from synthetic_neowise import EPOCH, field_file, target_file

from hd170312 import data, neowise, preregistered


def test_the_query_is_a_documented_tap_cone_of_ten_arcsec_with_the_columns_of_the_reading_and_sorted_by_time():
    url = data.neowise_target_url()
    assert url.startswith("https://irsa.ipac.caltech.edu/TAP/sync?") and "neowiser_p1bs_psd" in url and "CIRCLE" in url.upper().replace("+", " ") and "FORMAT=csv" in url
    for column in data.NEOWISE_COLUMNS.split(","):
        assert column in url
    assert "order+by+mjd" in url.lower() and "0.002777" in url


def test_the_cuts_drop_a_frame_of_quality_zero_a_flagged_band_a_poor_photometry_and_a_source_six_arcsec_away(tmp_path):
    rows = neowise.read_exposures(target_file(tmp_path / "t.csv"))
    assert rows["btjd"].size == 28 * 13, "the four rows that each fail one cut are dropped, the 364 others stay"
    assert np.all(np.diff(rows["btjd"]) > 0)


def test_the_visits_are_runs_with_no_gap_of_a_day_and_the_visit_with_the_window_holds_seven_of_its_thirteen_exposures(tmp_path):
    rows = neowise.read_exposures(target_file(tmp_path / "t.csv"))
    assert rows["visit"].max() + 1 == 28 and np.bincount(rows["visit"]).tolist() == [13] * 28
    window = preregistered.only_window(rows["btjd"])
    assert set(rows["visit"][window]) == {2} and window.sum() == 7, "the window holds exposures of one visit only"


def test_the_time_is_bjd_tdb_minus_2457000_within_seven_minutes_of_the_mjd_and_the_two_differ(tmp_path):
    rows = neowise.read_exposures(target_file(tmp_path / "t.csv"))
    offset = (rows["btjd"] - (rows["mjd"] + 2400000.5 - 2457000.0)) * 1440.0
    assert np.abs(offset).max() < 7.1 and np.abs(offset).max() > 1.0, "TDB - UTC is about a minute and the light-travel time to the barycentre is at most 5.8 minutes for this star (ecliptic latitude -45 deg)"


def test_the_numbers_of_the_visit_give_the_window_dip_the_minimum_time_and_its_distance_from_the_prediction(tmp_path):
    numbers = neowise.visit_numbers(neowise.read_exposures(target_file(tmp_path / "t.csv")))
    assert numbers["visits"] == 28 and numbers["window_exposures"] == 7 and numbers["other_exposures"] == 6
    assert 0.04 < numbers["w2_window_minus_other"] < 0.12 and abs(numbers["w2_minimum_time"] - EPOCH) < 0.25
    assert numbers["w2_minimum_offset_days"] == abs(numbers["w2_minimum_time"] - numbers["predicted_time"]) and abs(numbers["predicted_time"] - EPOCH) < 0.01


def test_a_flat_series_gives_a_window_difference_of_zero(tmp_path):
    numbers = neowise.visit_numbers(neowise.read_exposures(target_file(tmp_path / "t.csv", dip=0.0)))
    assert abs(numbers["w2_window_minus_other"]) < 0.03


def test_the_comparison_stars_are_the_bright_sources_with_a_counterpart_and_the_star_and_its_blend_are_not_among_them(tmp_path):
    rows = neowise.read_exposures(target_file(tmp_path / "t.csv"))
    shifts = neowise.comparison_shifts(rows, field_file(tmp_path / "f.csv", [0.0, 0.02, -0.01, 0.005]))
    assert sorted(round(v, 2) for v in shifts.values()) == [-0.01, 0.0, 0.01 - 0.005, 0.02] or len(shifts) == 4, shifts
    assert len(shifts) == 4 and all("1" != name and "98" != name and "99" != name and name != "" for name in shifts)
    assert abs(max(shifts.values()) - 0.02) < 0.01 and abs(min(shifts.values()) + 0.01) < 0.01


def test_a_comparison_star_that_dips_with_the_window_is_the_largest_shift(tmp_path):
    rows = neowise.read_exposures(target_file(tmp_path / "t.csv"))
    shifts = neowise.comparison_shifts(rows, field_file(tmp_path / "f.csv", [0.0, 0.0, 0.15]))
    assert max(shifts.values()) > 0.12


def test_the_dip_of_a_visit_is_its_faintest_run_of_eight_exposures_against_the_rest_and_a_visit_of_fewer_than_twelve_is_not_tested():
    w2 = np.full(12, 7.5)
    w2[3:11] += 0.1
    assert abs(neowise.window_dip(w2) - 0.1) < 1e-12, "eight exposures 0.1 mag fainter than the other four"
    assert abs(neowise.window_dip(np.full(12, 7.5))) < 1e-12 and neowise.window_dip(np.array([7.5] * 11 + [7.6])) > 0
    rows = {"visit": np.array([0] * 12 + [1] * 11), "w2": np.concatenate([w2, np.full(11, 7.5)]), "mjd": np.arange(23) * 0.0625}
    assert neowise.dips_of_visits(rows) == {0: neowise.window_dip(w2)}, "the visit of 11 exposures is left out"


def test_the_null_of_the_visits_has_the_eclipse_visit_far_above_the_others_and_gives_the_spacing_of_its_exposures(tmp_path):
    rows = neowise.read_exposures(target_file(tmp_path / "t.csv"))
    null = neowise.dip_null(rows)
    assert null["visits_tested"] == 28 and null["eclipse_visit_dip"] > 0.04 and null["largest_other_dip"] < 0.04 and null["eclipse_visit_dip"] > 1.5 * null["largest_other_dip"] and null["largest_other_dip"] >= null["median_other_dip"]
    assert abs(null["eclipse_visit_spacing_minutes"] - 144.0) < 1.0, "the synthetic exposures are 0.1 d apart"
    assert 0.01 < null["median_visit_sd"] < 0.03
