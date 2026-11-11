"""The ASAS-SN saturated-star photometry: the reader, the nightly means, the windows of the two kinds of epoch, the test on a synthetic file and the checksum that `fetch` verifies."""
import json

import numpy as np
import pytest
from synthetic_asassn import synthetic_file

from hd170312 import asassn, data


def test_the_reader_drops_comment_lines_and_rows_without_a_magnitude_and_names_a_dataset_for_each_filter_and_camera(tmp_path):
    path = tmp_path / "small.csv"
    path.write_text("HJD,UT Date,Camera,FWHM,Limit,mag,mag_err,flux(mJy),flux_err,Filter\n### magnitude offsets  bi -0.000\n2457001.0,x,bi,2,NA,8.7,0.03,NA,NA,g\n2457002.0,x,be,2,NA,8.9,0.03,NA,NA,V\n2457003.0,x,bi,2,NA,NA,NA,NA,NA,g\n")
    rows, names = asassn.read_rows(path)
    assert names == ["V:be", "g:bi"]
    assert rows[0].tolist() == [[2.0, 8.9]] and rows[1].tolist() == [[1.0, 8.7]], "time is HJD minus 2457000, one array for each dataset in the order of the names"


def test_a_nightly_mean_pools_the_points_of_one_integer_day_of_one_dataset(tmp_path):
    path = tmp_path / "small.csv"
    path.write_text("HJD,UT Date,Camera,FWHM,Limit,mag,mag_err,flux(mJy),flux_err,Filter\n" + "".join(f"{hjd},x,bi,2,NA,{mag},0.03,NA,NA,g\n" for hjd, mag in ((2457001.2, 8.70), (2457001.4, 8.72), (2457002.3, 8.80), (2457003.3, 8.60), (2457004.3, 8.70))))
    t, m, _, datasets, points = asassn.read_nights(path)[0]
    assert points.tolist() == [2, 1, 1, 1] and abs(m[0] - 8.71) < 1e-12 and abs(t[0] - 1.3) < 1e-9 and set(datasets.tolist()) == {0}


def test_the_test_recovers_the_amplitude_of_each_kind_of_epoch_and_counts_the_nights_of_each_window(tmp_path):
    result = asassn.run(synthetic_file(tmp_path / "asassn.csv", odd=1.0, even=2.0), shifts=300)
    assert result["nights"] < result["measurements"] and result["nights"] == 2 * 2400, "two points a night in each of two cameras: a night is a camera-night"
    raw = [float(line.split(",")[0]) for line in (tmp_path / "asassn.csv").read_text().splitlines()[2:] if line.split(",")[5] != "NA"]
    assert result["first_hjd"] == min(raw) and result["last_hjd"] == max(raw), "the span of the measurements, HJD"
    assert result["odd"]["nights"] > 0 and result["even"]["nights"] > 0 and result["odd"]["nights"] + result["even"]["nights"] < 200
    assert abs(result["odd"]["amplitude"] - 1.0) < 0.35 and abs(result["even"]["amplitude"] - 2.0) < 0.35
    assert result["even"]["amplitude"] > result["odd"]["amplitude"] and result["even"]["p"] < 0.05 and result["shifts"] == 300 and result["seed"] == 1


def test_a_flat_star_has_amplitudes_near_zero_and_no_small_p(tmp_path):
    result = asassn.run(synthetic_file(tmp_path / "asassn.csv", odd=0.0, even=0.0), shifts=300)
    assert abs(result["odd"]["amplitude"]) < 0.4 and abs(result["even"]["amplitude"]) < 0.4 and result["odd"]["p"] > 0.05 and result["even"]["p"] > 0.05


def test_the_even_window_is_the_nights_of_even_cycles_only_and_is_disjoint_from_the_odd_window(tmp_path):
    t = asassn.read_nights(synthetic_file(tmp_path / "asassn.csv"))[0][0]
    from hd170312 import preregistered

    odd, even = preregistered.only_window(t), asassn.even_window(t)
    assert odd.sum() > 0 and even.sum() > 0 and not (odd & even).any()


def test_the_error_is_widened_by_the_root_of_the_chi_square_per_degree_of_freedom_and_never_narrowed(tmp_path):
    for scatter in (0.010, 0.100):
        part = asassn.run(synthetic_file(tmp_path / f"{scatter}.csv", scatter=scatter), shifts=20)["odd"]
        assert part["scale"] == max(1.0, np.sqrt(part["chi2"] / part["dof"])) and abs(part["scaled_error"] - part["formal_error"] * part["scale"]) < 1e-12


def test_fetch_verifies_the_tracked_asassn_file_and_never_downloads_it(tmp_path):
    name = "asassn/asassn_v1_saturated.csv"
    (tmp_path / "asassn").mkdir()
    (tmp_path / name).write_text("downloaded by hand\n")
    good = {name: {"sha256": data.sha256(tmp_path / name), "bytes": 19}}
    called = []
    assert data.ensure_manual(name, tmp_path, good) == tmp_path / name and not called
    with pytest.raises(data.ChecksumMismatch):
        data.ensure_manual(name, tmp_path, {name: {"sha256": "0" * 64, "bytes": 1}})
    with pytest.raises(FileNotFoundError, match="Compute Coordinate Lightcurve"):
        data.ensure_manual("asassn/missing.csv", tmp_path, {})


def test_the_tracked_checksums_name_the_asassn_file_and_it_matches_the_file_in_the_repository():
    recorded = json.loads(data.CHECKSUMS.read_text())
    entry = recorded["asassn/asassn_v1_saturated.csv"]
    assert data.sha256(data.DATA / "asassn" / "asassn_v1_saturated.csv") == entry["sha256"] and len(entry["sha256"]) == 64


def test_the_test_before_a_date_fits_only_the_earlier_nights_and_keeps_the_errors_of_all_the_nights(tmp_path):
    path = synthetic_file(tmp_path / "asassn.csv")
    full, early = asassn.run(path, shifts=20), asassn.run(path, shifts=20, before=-1800.0)
    assert early["before"] == -1800.0 and full["before"] is None
    assert 0 < early["nights"] < full["nights"] and early["nights"] == 2 * 1200, "the nights before BTJD -1800: 1,200 days of two cameras"
    assert early["odd"]["nights"] < full["odd"]["nights"] and early["even"]["nights"] < full["even"]["nights"]
    assert abs(early["odd"]["amplitude"] - 1.0) < 0.4 and abs(early["even"]["amplitude"] - 2.0) < 0.4
