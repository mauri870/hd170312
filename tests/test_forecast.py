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
