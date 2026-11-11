"""The forecast of the deeper eclipse in sector 114."""
from hd170312 import forecast


def test_the_forecast_is_the_deeper_eclipse_of_the_shortest_period_inside_sector_114_and_no_other_period_has_one():
    events = [{"time": 2377.1543, "error": 0.00015, "kind": "shallow"}, {"time": 4131.0663, "error": 0.00012, "kind": "deep"}, {"time": 4169.8482, "error": 0.0005, "kind": "shallow"}]
    rows = forecast.forecast(events)
    deciding = forecast.deciding_eclipse(rows)
    assert deciding["N"] == 26 and abs(deciding["time"] - 4475.8152) < 2e-4 and deciding["error"] < 3e-4 and deciding["date"] == "2027-03-11"
    assert abs(deciding["days_before_the_end"] - 2.68) < 0.01
    assert all(times == [] for times in deciding["deeper_eclipses_in_the_window_of_the_other_periods"].values())
    assert forecast.minimum_gap(rows)[26] == 1 and sum(forecast.minimum_gap(rows).values()) == 1


EVENTS = [{"time": 2377.1543, "error": 0.00015, "kind": "shallow"}, {"time": 4131.0663, "error": 0.00012, "kind": "deep"}, {"time": 4169.8482, "error": 0.0005, "kind": "shallow"}]


def test_the_eccentricity_bound_comes_from_the_phase_of_the_deep_eclipse_at_each_period():
    bounds = forecast.eccentricity(forecast.forecast(EVENTS))
    assert abs(bounds[26]["phase_of_deep"] - 0.4375) < 1e-3 and abs(bounds[26]["e_cos_omega_at_least"] - 0.098) < 2e-3
    assert abs(bounds[13]["phase_of_deep"] - 0.7188) < 1e-3 and abs(bounds[13]["e_cos_omega_at_least"] - 0.344) < 2e-3


def test_a_secondary_at_half_phase_gives_no_eccentricity_bound():
    events = [{"time": 0.0, "error": 1e-4, "kind": "shallow"}, {"time": 5.0, "error": 1e-4, "kind": "deep"}, {"time": 100.0, "error": 1e-4, "kind": "shallow"}]
    assert forecast.eccentricity(forecast.forecast(events, divisors=(10,)))[10]["e_cos_omega_at_least"] < 1e-9


def test_the_eclipses_that_only_68_95_d_predicts_before_sector_114_include_2026_december_2():
    rows = forecast.only_shortest(EVENTS, low=forecast.btjd_of("2026-11-11"), high=forecast.SECTOR_114[0])
    shallow = [row for row in rows if row["kind"] == "shallow"]
    assert shallow and shallow[0]["date"] == "2026-12-02" and abs(shallow[0]["time"] - 4376.6975) < 2e-3
    assert all(row["kind"] in ("shallow", "deep") and not row["also_137_90"] for row in rows)
