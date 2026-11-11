"""The forecast of the deeper eclipse in TESS sector 114.

For each period P = DeltaT / N the three complete eclipses give an ephemeris, and the eclipse of each kind in the window of sector 114 is predicted with its error. The window, BTJD 4449.5 to 4478.5, is the planned
dates of the sector (https://heasarc.gsfc.nasa.gov/docs/tess/sector.html)."""
import numpy as np

from hd170312.eclipses import assign_cycles, fit_ephemeris, predict

SECTOR_114 = (4449.5, 4478.5)  # BTJD at 0 h UT of the planned start and stop
DIVISORS = (1, 2, 3, 4, 9, 13, 26)
JD_OFFSET = 2457000.0


def forecast(events, divisors=DIVISORS, window=SECTOR_114):
    """{N: {period, period_error, epochs, predictions}} for P = (T_3 - T_1) / N, with the predicted eclipses within 6 d before the window and 1 d after it."""
    shallow = [e for e in events if e["kind"] == "shallow"]
    separation = shallow[-1]["time"] - shallow[0]["time"]
    rows = {}
    for divisor in divisors:
        period = separation / divisor
        fit_events = [{"time": e["time"], "error": e["error"], "kind": e["kind"]} for e in events]
        for event, cycle in zip(fit_events, assign_cycles([e["time"] for e in fit_events], [e["kind"] for e in fit_events], period)):
            event["cycle"] = cycle
        fit = fit_ephemeris(fit_events)
        predictions = []
        for kind in fit["kinds"]:
            for cycle in range(-1, 40):
                time, error = predict(fit, kind, cycle)
                if window[0] - 6 <= time <= window[1] + 1:
                    predictions.append({"kind": kind, "cycle": cycle, "time": time, "error": error})
        rows[divisor] = {"period": fit["period"], "period_error": fit["period_error"], "epochs": fit["epochs"], "predictions": predictions}
    return rows


def calendar_date(btjd):
    """The calendar date (UTC) of a BTJD, from BJD = BTJD + 2457000 read as a Julian date of the TDB scale."""
    from astropy.time import Time

    return Time(btjd + JD_OFFSET, format="jd", scale="tdb").utc.strftime("%Y-%m-%d")


def deciding_eclipse(rows, window=SECTOR_114):
    """The deeper-kind eclipse that only the shortest period predicts in the window: (N, prediction) of the largest N, with the days before the end of the window and the JD and date."""
    divisor = max(rows)
    inside = [p for p in rows[divisor]["predictions"] if p["kind"] == "deep" and window[0] <= p["time"] <= window[1]]
    if not inside:
        return None
    prediction = inside[0]
    others = {n: [p["time"] for p in row["predictions"] if p["kind"] == "deep" and window[0] <= p["time"] <= window[1]] for n, row in rows.items() if n != divisor}
    return {"N": divisor, **prediction, "jd": prediction["time"] + JD_OFFSET, "date": calendar_date(prediction["time"]), "days_before_the_end": float(window[1] - prediction["time"]), "deeper_eclipses_in_the_window_of_the_other_periods": {n: [float(t) for t in times] for n, times in others.items()}}


def minimum_gap(rows, window=SECTOR_114):
    """For each period, the number of its predicted eclipses (either kind) inside the window."""
    return {n: int(np.sum([window[0] <= p["time"] <= window[1] for p in row["predictions"]])) for n, row in rows.items()}
