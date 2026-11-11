"""The ASAS-3 light curve against the eclipse trains of the seven periods: the pooled matched filter with its period-preserving shift null.

The light curve is grades A and B in the aperture of lowest median frame error, as relative flux by dataset. For each period the model depth at a point is the TESS trapezoid of each kind at each predicted epoch. The
amplitude A is the sum of m y over the sum of m squared (m the model, y the deficit against the median within 15 d); A = 1 is the TESS depth. The null shifts the train by one offset between 1 d and P - 1 d."""
import math

import numpy as np

from hd170312.forecast import forecast

SHAPES = {"shallow": (0.0855, 0.747, 0.785), "deep": (0.1536, 0.653, 0.811)}  # flux depth, duration (d), ingress share of the three complete eclipses
SEED = 952
SHIFTS = 1000
REFERENCE_DAYS = 15.0
NS = (1, 2, 3, 4, 9, 13, 26)
APERTURES = ("MAG_0", "MAG_1", "MAG_2", "MAG_3", "MAG_4")


def prepare(rows):
    """(time in BTJD, relative flux, aperture, rows read) of the points of grades A and B in the aperture of the lowest median frame error."""
    aperture = min(APERTURES, key=lambda a: np.median([r["err"][a] for r in rows]))
    pairs = []
    for dataset in sorted({r["dataset"] for r in rows}):
        mine = [r for r in rows if r["dataset"] == dataset and r["grade"] in "AB" and r["mag"][aperture] < 29.0]
        mags = np.array([r["mag"][aperture] for r in mine])
        pairs += [(r["hjd"] - 7000.0, 10 ** (-0.4 * (m - np.median(mags)))) for r, m in zip(mine, mags, strict=True)]
    pairs.sort()
    return np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs]), aperture, len(rows)


def ephemerides(events):
    """{N: (period, {kind: epoch})}: the ephemeris of each period fitted to the three complete eclipses [{time, error, kind}] (hd170312.forecast)."""
    return {n: (row["period"], row["epochs"]) for n, row in forecast(events, NS).items()}


def predicted(eph, lo, hi):
    """{N: [(kind, time)]}: the predicted eclipses of each period within [lo, hi] (BTJD)."""
    return {n: [(kind, epoch + k * period) for kind, epoch in epochs.items() for k in range(math.ceil((lo - epoch) / period), math.floor((hi - epoch) / period) + 1)] for n, (period, epochs) in eph.items()}


def model(time, events):
    """The model depth at each point: the trapezoid of each event's kind, the largest where events overlap."""
    m = np.zeros(time.size)
    for kind, t in events:
        depth, duration, ingress = SHAPES[kind]
        half = duration / 2
        m = np.maximum(m, np.clip((half - np.abs(time - t)) / (ingress * half), 0.0, 1.0) * depth)
    return m


def reference(time, flux, excluded):
    """The median flux within REFERENCE_DAYS of each point of the points farther than 1 d from every predicted eclipse."""
    keep = np.all(np.abs(time[:, None] - excluded[None, :]) > 1.0, axis=1)
    out = np.empty(time.size)
    for i, t in enumerate(time):
        near = keep & (np.abs(time - t) <= REFERENCE_DAYS)
        out[i] = np.median(flux[near]) if near.sum() >= 3 else np.median(flux[keep])
    return out


def amplitude(time, deficit, events, scatter):
    """(A, its error, points in the model) or None when no point lies in the model."""
    m = model(time, events)
    used = m > 0
    if not used.any():
        return None
    s2 = float(np.sum(m[used] ** 2))
    return float(np.sum(m[used] * deficit[used]) / s2), scatter / math.sqrt(s2), int(used.sum())


class Light:
    """The light curve with the deficits against the trains of all periods and the scatter of the points that no model covers."""

    def __init__(self, time, flux, eph):
        self.time, self.flux, self.eph = time, flux, eph
        self.events = predicted(eph, float(time.min()) - 1.0, float(time.max()) + 1.0)
        self.excluded = np.array(sorted({t for e in self.events.values() for _, t in e}))
        self.deficit = 1.0 - flux / reference(time, flux, self.excluded)
        self.none = np.all(np.abs(time[:, None] - self.excluded[None, :]) > 1.0, axis=1)  # points with no model at any period
        self.scatter = float(np.std(self.deficit[self.none]))


def shifted_amplitudes(light, selection, period, rng, shifts):
    """The amplitudes (and the amplitude over its error) of `shifts` shifted copies of a train: one offset between 1 d and P - 1 d for the whole train."""
    amplitudes, ratios = np.zeros(shifts), np.zeros(shifts)
    for j in range(shifts):
        offset = float(rng.uniform(1.0, period - 1.0))
        found = amplitude(light.time, light.deficit, [(k, t + offset) for k, t in selection], light.scatter)
        if found:
            amplitudes[j], ratios[j] = found[0], found[0] / found[1]
    return amplitudes, ratios


def shift_test(light, selection, period, rng, shifts):
    """{events, points, A, sigma_A, p_A, p_A_over_sigma} of a train against `shifts` shifted trains, p = (1 + shifted at least as large) / (1 + shifts)."""
    found = amplitude(light.time, light.deficit, selection, light.scatter)
    if found is None:
        return None
    a0, sigma, points = found
    amplitudes, ratios = shifted_amplitudes(light, selection, period, rng, shifts)
    return {"events": len(selection), "points": points, "A": a0, "sigma_A": sigma, "p_A": float((1 + np.sum(amplitudes >= a0)) / (1 + shifts)), "p_A_over_sigma": float((1 + np.sum(ratios >= a0 / sigma)) / (1 + shifts))}


def subset_tests(light, shifts=50_000, seed=20261109):
    """The shift tests at `shifts` shifts of the train of 68.95 d (all events), of its events that 137.90 d does not predict, of 137.90 d, and of the events that only 199.19 d predicts; the generator is drawn in this order."""
    rng = np.random.default_rng(seed)
    events = light.events
    others13 = np.array([t for _, t in events[13]])
    unique26 = [(k, t) for k, t in events[26] if not np.any(np.abs(others13 - t) <= 0.5)]
    others_than9 = np.array([t for n, e in events.items() if n != 9 for _, t in e])
    unique9 = [(k, t) for k, t in events[9] if not np.any(np.abs(others_than9 - t) <= 0.5)]
    return {"N26_all": shift_test(light, events[26], light.eph[26][0], rng, shifts), "N26_events_that_N13_does_not_predict": shift_test(light, unique26, light.eph[26][0], rng, shifts),
            "N13_all": shift_test(light, events[13], light.eph[13][0], rng, shifts), "N9_events_no_other_period_predicts": shift_test(light, unique9, light.eph[9][0], rng, shifts), "shifts": shifts}
