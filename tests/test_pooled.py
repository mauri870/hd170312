"""The ASAS-3 pooled matched filter on a synthetic archive file."""
import numpy as np
from synthetic_asas import synthetic_file

from hd170312 import asas, pooled

EVENTS = [{"time": 2377.1543, "error": 0.00015, "kind": "shallow"}, {"time": 4131.0663, "error": 0.00012, "kind": "deep"}, {"time": 4169.8482, "error": 0.0005, "kind": "shallow"}]


def light_of(path):
    time, flux, aperture, rows = pooled.prepare(asas.read_table(path))
    return pooled.Light(time, flux, pooled.ephemerides(EVENTS)), aperture, rows


def test_columns_are_read_by_their_names_when_the_datasets_list_the_apertures_in_another_order(tmp_path):
    text = """# HJD MAG_3 MAG_0 MAG_1 MAG_2 MAG_4 MER_3 MER_0 MER_1 MER_2 MER_4 GRADE FRAME
#dataset= 1 ; 1 F1
#desig= 183419-6849.9
 1870.5 8.30 8.10 8.20 8.40 8.50 0.01 0.05 0.04 0.03 0.02  A 1
# HJD MAG_4 MAG_0 MAG_1 MAG_2 MAG_3 MER_4 MER_0 MER_1 MER_2 MER_3 GRADE FRAME
#dataset= 2 ; 1 F2
#desig= 183419-6849.9
 1872.5 8.50 8.10 8.20 8.40 8.30 0.02 0.05 0.04 0.03 0.01  B 2
"""
    path = tmp_path / "a.txt"
    path.write_text(text)
    rows = asas.read_table(path)
    assert [r["mag"]["MAG_3"] for r in rows] == [8.30, 8.30] and [r["mag"]["MAG_4"] for r in rows] == [8.50, 8.50]
    time, _flux, aperture, count = pooled.prepare(rows)
    assert aperture == "MAG_3" and count == 2 and time.tolist() == [-5129.5, -5127.5]


def test_the_pooled_amplitude_of_68_95_d_is_found_on_a_train_injected_at_its_epochs_and_not_on_a_flat_curve(tmp_path):
    injected, _, _ = light_of(synthetic_file(tmp_path / "injected.txt", amplitude=1.0))
    flat, _, _ = light_of(synthetic_file(tmp_path / "flat.txt", amplitude=0.0))
    a_injected = pooled.amplitude(injected.time, injected.deficit, injected.events[26], injected.scatter)
    a_flat = pooled.amplitude(flat.time, flat.deficit, flat.events[26], flat.scatter)
    assert a_injected[0] > 0.5 and a_injected[0] / a_injected[1] > 4 and abs(a_flat[0]) < 4 * a_flat[1]


def test_the_shift_test_counts_its_probability_as_written_and_is_reproducible(tmp_path):
    light, _, _ = light_of(synthetic_file(tmp_path / "a.txt"))
    first = pooled.shift_test(light, light.events[26], light.eph[26][0], np.random.default_rng(5), 60)
    second = pooled.shift_test(light, light.events[26], light.eph[26][0], np.random.default_rng(5), 60)
    assert first == second and first["p_A"] * 61 == round(first["p_A"] * 61) and first["p_A"] >= 1 / 61
