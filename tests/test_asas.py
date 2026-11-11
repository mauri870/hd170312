"""The ASAS-3 reading and the pre-registered test, on a synthetic archive file."""
import numpy as np
import pytest
from synthetic_asas import synthetic_file

from hd170312 import asas, preregistered


def test_the_reader_keeps_grade_a_measured_rows_of_the_designation_and_merges_blocks_by_dataset(tmp_path):
    text = """#dataset= 1 ; 1 F1
#desig= 183419-6849.9
# HJD MAG_3 MAG_0 MAG_1 MAG_2 MAG_4 MER_3 MER_0 MER_1 MER_2 MER_4 GRADE FRAME
 1870.50 8.700 8.7 8.7 8.7 8.7 0.03 0.03 0.03 0.03 0.03  A 1
 1872.50 8.710 8.7 8.7 8.7 8.7 0.03 0.03 0.03 0.03 0.03  B 2
 1874.50 29.999 8.7 8.7 8.7 8.7 0.03 0.03 0.03 0.03 0.03  A 3
#dataset= 1 ; 2 F1
#desig= 183416-6849.6
 1876.50 9.000 8.7 8.7 8.7 8.7 0.03 0.03 0.03 0.03 0.03  A 4
#dataset= 2 ; 1 F2
#desig= 183419-6849.9
# HJD MAG_3 MAG_0 MAG_1 MAG_2 MAG_4 MER_3 MER_0 MER_1 MER_2 MER_4 GRADE FRAME
 1880.50 8.720 8.7 8.7 8.7 8.7 0.03 0.03 0.03 0.03 0.03  A 5
"""
    path = tmp_path / "a.txt"
    path.write_text(text)
    rows = asas.read_rows(path)
    assert sorted(rows) == [1, 2] and rows[1].tolist() == [[-5129.5, 8.7]] and rows[2].tolist() == [[-5119.5, 8.72]]


def test_nights_are_the_means_of_the_points_of_a_day_with_the_scatter_of_the_means_over_the_root_of_the_points():
    rows = {1: np.array([[-5000.2, 8.0], [-5000.4, 8.2], [-4998.5, 8.1], [-4996.5, 8.3], [-4994.5, 8.1]])}
    nights = asas.nightly_means(rows)[1]
    assert list(nights["c"]) == [2, 1, 1, 1] and nights["m"][0] == pytest.approx(8.1) and nights["t"][0] == pytest.approx(-5000.3)
    assert nights["e"][0] == nights["sc"] / np.sqrt(2) and nights["e"][1] == nights["sc"]


def test_an_injected_eclipse_train_is_recovered_by_the_test_and_a_flat_light_curve_is_not_preferred(tmp_path):
    injected = preregistered.run(synthetic_file(tmp_path / "injected.txt", amplitude=1.0), shifts=300)
    assert injected["nights"] >= 6 and abs(injected["amplitude"] - 1.0) < 4 * injected["formal_error"] + 0.15
    assert injected["decision"] == "favours 68.95 d"
    flat = preregistered.run(synthetic_file(tmp_path / "flat.txt", amplitude=0.0), shifts=300)
    assert abs(flat["amplitude"]) < 4 * flat["formal_error"] and flat["decision"] != "favours 68.95 d"


def test_the_null_is_reproducible_with_its_seed_and_its_p_is_counted_as_written(tmp_path):
    arrays = asas.stacked(asas.nightly_means(asas.read_rows(synthetic_file(tmp_path / "a.txt"))))
    first = preregistered.null_distribution(arrays, 50, seed=1)
    assert np.array_equal(first, preregistered.null_distribution(arrays, 50, seed=1)) and not np.array_equal(first, preregistered.null_distribution(arrays, 50, seed=2))
    result = preregistered.run(tmp_path / "a.txt", shifts=50, seed=1)
    assert result["p"] == (1 + result["shifts_at_least_as_large"]) / (1 + result["shifts"])



def test_a_null_with_no_exclusion_keeps_every_shift_the_generator_draws_and_the_run_reports_its_p(tmp_path):
    arrays = asas.stacked(asas.nightly_means(asas.read_rows(synthetic_file(tmp_path / "a.txt"))))
    unrestricted = preregistered.null_distribution(arrays, 40, seed=3, exclusion=0.0)
    draws = np.random.default_rng(3).uniform(0, preregistered.P26, 40)
    assert np.allclose(unrestricted, [preregistered.a_only(arrays, u)[0][0] for u in draws]), "no shift is rejected: the first 40 draws are the null"
    result = preregistered.run(tmp_path / "a.txt", shifts=40, seed=3)
    assert result["p_no_exclusion"] == (1 + int((unrestricted >= result["amplitude"]).sum())) / 41 and "null_amplitudes_no_exclusion" not in result


def test_the_rule_fixed_in_advance_has_three_outcomes():
    assert preregistered.decision(0.0066, 1.035) == "favours 68.95 d"
    assert preregistered.decision(0.005, 1.8) == "weak preference"
    assert preregistered.decision(0.2, 1.0) == "no preference"
    assert preregistered.decision(0.05, 1.0) == "weak preference"



def test_read_rows_takes_the_aperture_by_its_name_when_the_datasets_list_the_columns_in_another_order(tmp_path):
    """The archive lists MAG_3 first in some datasets of a star and MAG_4 first in others; a reader that takes the second column reads two different apertures."""
    text = """#dataset= 1 ; 1 F1
#desig= 183419-6849.9
#     HJD      MAG_3  MAG_0  MAG_1  MAG_2  MAG_4    MER_3 MER_0 MER_1 MER_2 MER_4 GRADE FRAME
 1870.50  8.300  8.100  8.200  8.400  8.500  0.01 0.05 0.04 0.03 0.02  A 1
#dataset= 2 ; 1 F2
#desig= 183419-6849.9
#     HJD      MAG_4  MAG_0  MAG_1  MAG_2  MAG_3    MER_4 MER_0 MER_1 MER_2 MER_3 GRADE FRAME
 1872.50  8.500  8.100  8.200  8.400  8.300  0.02 0.05 0.04 0.03 0.01  A 2
"""
    path = tmp_path / "a.txt"
    path.write_text(text)
    assert asas.read_rows(path)[1].tolist() == [[-5129.5, 8.3]] and asas.read_rows(path)[2].tolist() == [[-5127.5, 8.3]]
    assert asas.read_rows(path, aperture="MAG_4")[2].tolist() == [[-5127.5, 8.5]]


def test_the_result_of_the_test_keeps_the_null_amplitudes_for_the_figure(tmp_path):
    found = preregistered.run(synthetic_file(tmp_path / "a.txt"), shifts=40)
    assert len(found["null_amplitudes"]) == 40 == found["shifts"] and all(isinstance(value, float) for value in found["null_amplitudes"])
