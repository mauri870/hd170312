"""The forecast of the deeper eclipse in TESS sector 114.

For each period P = DeltaT / N the three complete eclipses give an ephemeris, and the eclipse of each kind in the window of sector 114 is predicted with its error. The window, BTJD 4449.5 to 4478.5, is the planned
dates of the sector (https://heasarc.gsfc.nasa.gov/docs/tess/sector.html)."""
import numpy as np

from hd170312.eclipses import assign_cycles, fit_ephemeris, predict

SECTOR_114 = (4449.5, 4478.5)  # BTJD at 0 h UT of the planned start and stop
DIVISORS = (1, 2, 3, 4, 9, 13, 26)
JD_OFFSET = 2457000.0
NOTE_DATE = "2026-11-11"  # the date of the Note; the eclipses that only 68.95 d predicts are listed from it to sector 114


def _fit(events, divisor):
    """The ephemeris fitted to the events at P = (T_3 - T_1) / divisor, the separation of the first and last shallow eclipses."""
    shallow = [e for e in events if e["kind"] == "shallow"]
    period = (shallow[-1]["time"] - shallow[0]["time"]) / divisor
    fit_events = [{"time": e["time"], "error": e["error"], "kind": e["kind"]} for e in events]
    for event, cycle in zip(fit_events, assign_cycles([e["time"] for e in fit_events], [e["kind"] for e in fit_events], period)):
        event["cycle"] = cycle
    return fit_ephemeris(fit_events)


def forecast(events, divisors=DIVISORS, window=SECTOR_114):
    """{N: {period, period_error, epochs, predictions}} for P = (T_3 - T_1) / N, with the predicted eclipses within 6 d before the window and 1 d after it."""
    rows = {}
    for divisor in divisors:
        fit = _fit(events, divisor)
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


def btjd_of(date):
    """The BTJD of 0 h UTC of a calendar date, the inverse of calendar_date."""
    from astropy.time import Time

    return float(Time(date, scale="utc").tdb.jd - JD_OFFSET)


def eccentricity(rows):
    """{N: {phase_of_deep, e_cos_omega_at_least}}: the phase of the deep eclipse after the shallow one in each period's fit, and the lower bound on the eccentricity it implies, |e cos w| = (pi / 2) |phase - 1/2|, to first order in e (Winn 2010, in Exoplanets, ed. Seager, eq. 33: the secondary eclipse falls at P / 2 (1 + 4 e cos w / pi) after the primary)."""
    bounds = {}
    for divisor, row in rows.items():
        phase = ((row["epochs"]["deep"] - row["epochs"]["shallow"]) / row["period"]) % 1.0
        bounds[divisor] = {"phase_of_deep": float(phase), "e_cos_omega_at_least": float(np.pi / 2 * abs(phase - 0.5))}
    return bounds


def only_shortest(events, low, high, shortest=26, longer=13):
    """The eclipses of each kind that P = separation / shortest predicts between BTJD low and high and P = separation / longer does not (none of its eclipses of the same kind within 0.5 d), with their errors and UTC dates and times."""
    from astropy.time import Time

    short, long_ = _fit(events, shortest), _fit(events, longer)
    rows = []
    for kind in short["kinds"]:
        others = [long_["epochs"][kind] + cycle * long_["period"] for cycle in range(-200, 200)]
        first = int(np.ceil((low - short["epochs"][kind]) / short["period"]))
        last = int(np.floor((high - short["epochs"][kind]) / short["period"]))
        for cycle in range(first, last + 1):
            time, error = predict(short, kind, cycle)
            also = any(abs(time - other) < 0.5 for other in others)
            if not also:
                utc = Time(time + JD_OFFSET, format="jd", scale="tdb").utc
                rows.append({"kind": kind, "cycle": cycle, "time": time, "error": error, "date": utc.strftime("%Y-%m-%d"), "utc": utc.strftime("%Y-%m-%d %H:%M"), "also_137_90": also})
    return sorted(rows, key=lambda row: row["time"])


def minimum_gap(rows, window=SECTOR_114):
    """For each period, the number of its predicted eclipses (either kind) inside the window."""
    return {n: int(np.sum([window[0] <= p["time"] <= window[1] for p in row["predictions"]])) for n, row in rows.items()}
