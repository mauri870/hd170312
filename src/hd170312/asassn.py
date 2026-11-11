"""The ASAS-SN photometry of saturated stars (Winecki & Kochanek 2024; Shappee et al. 2014, ApJ 788, 48; Kochanek et al. 2017, PASP 129, 104502) of HD 170312: the file, its nightly means and the two tests on its nights.

The file is the download of the Sky Patrol v1 form with the method "Saturated Stars (Machine Learning)" (see `data.ensure_manual`), one row for each measurement with the columns HJD, Camera, mag, mag_err and Filter; a row with no magnitude is dropped. A dataset is a filter and a camera, named
"filter:camera". A night is a camera-night: the nightly mean of the measurements of one dataset on one integer HJD day, as `asas.nightly_means` makes it, whose error is the scatter of the dataset's nightly means over the root of the points of the night.
The tests are those of `preregistered` on the nights within half a duration of an epoch of 68.95 d of one parity and of none of the other parity: the odd window holds the epochs that only 68.95 d predicts, and the even window, the control, those that both
periods predict. The error of an amplitude is widened by the root of the chi-square per degree of freedom of its fit when that exceeds one."""
import csv

import numpy as np

from hd170312 import preregistered
from hd170312.asas import nightly_means, stacked
from hd170312.data import ASASSN_FILE, DATA, verified

JD_OFFSET = 2457000.0  # the file's HJD to the offset of BTJD


def read_rows(path=None):
    """({dataset index: array of (HJD - 2457000, magnitude)}, dataset names): the rows with a finite magnitude and error; the names, 'filter:camera', are sorted and the index is the position in them."""
    path = verified(path) if path else verified(DATA / ASASSN_FILE)
    with open(path, newline="") as handle:
        table = [row for row in csv.DictReader(line for line in handle if not line.startswith("#"))]
    points = {}
    for row in table:
        if row["mag"] in ("NA", "") or row["mag_err"] in ("NA", ""):
            continue
        points.setdefault(f"{row['Filter']}:{row['Camera']}", []).append((float(row["HJD"]) - JD_OFFSET, float(row["mag"])))
    names = sorted(points)
    return {index: np.array(points[name]) for index, name in enumerate(names)}, names


def read_nights(path=None):
    """(the nights as the arrays of `asas.stacked` (time, magnitude, error, dataset index, points), the dataset names)."""
    rows, names = read_rows(path)
    return stacked(nightly_means(rows)), names


def even_window(t, shift=0.0):
    """The nights within half a duration of an even-cycle epoch and of no odd-cycle epoch: the epochs that both periods predict."""
    odd, even = np.zeros(len(t), bool), np.zeros(len(t), bool)
    for centre, _, duration, _ in preregistered.KINDS:
        k = np.round((t - centre - shift) / preregistered.P26)
        x = abs(t - (centre + shift + k * preregistered.P26))
        odd |= (x < duration / 2) & (k % 2 != 0)
        even |= (x < duration / 2) & (k % 2 == 0)
    return even & ~odd


def a_even(arrays, shift=0.0):
    """(A, its formal error, chi-square, offsets) of the even window, as `preregistered.a_only` does for the odd one."""
    t, m, e, datasets, _ = arrays
    window = even_window(t, shift)
    return preregistered.fit(t, m, e, datasets, np.where(window, preregistered.build(t, shift, "even"), 0.0)), window


def null_even(arrays, shifts=preregistered.SHIFTS, seed=preregistered.SEED):
    """The amplitudes of the even window at `shifts` shifted ephemerides, drawn as `preregistered.null_distribution` draws them."""
    generator = np.random.default_rng(seed)
    separation = (preregistered.KINDS[1][0] - preregistered.KINDS[0][0]) % preregistered.P26
    centres = (0.0, separation, (-separation) % preregistered.P26)
    values = []
    while len(values) < shifts:
        u = generator.uniform(0, preregistered.P26)
        if any(min((u - c) % preregistered.P26, (c - u) % preregistered.P26) < preregistered.EXCLUSION for c in centres):
            continue
        values.append(a_even(arrays, u)[0][0])
    return np.array(values)


def summary(fitted, window, null, nights):
    """The numbers of one window: its nights, A with its formal error and the error widened by the chi-square, the shifts that reach A and p = (1 + those) / (1 + shifts)."""
    (amplitude, sigma, chi2, _), dof = fitted, nights - (len(fitted[3]) + 1)
    scale = max(1.0, float(np.sqrt(chi2 / dof)))
    count = int((null >= amplitude).sum())
    return {"nights": int(window.sum()), "amplitude": float(amplitude), "formal_error": float(sigma), "chi2": float(chi2), "dof": int(dof), "scale": scale, "scaled_error": float(sigma * scale),
            "shifts_at_least_as_large": count, "p": (1 + count) / (1 + len(null)), "null_mean": float(null.mean()), "null_sd": float(null.std())}


def run(path=None, shifts=preregistered.SHIFTS, seed=preregistered.SEED, before=None):
    """The two tests on the file: the nights of all datasets, the odd window (the epochs that only 68.95 d predicts) and the even window (the control). With `before` (BTJD) only the nights earlier than that are fitted; the scatter that sets their errors is still that of
    all the nights of each dataset."""
    rows, names = read_rows(path)
    arrays = stacked(nightly_means(rows))
    if before is not None:
        arrays = tuple(column[arrays[0] < before] for column in arrays)
    nights = len(arrays[0])
    fitted, window = preregistered.a_only(arrays)
    odd = summary(fitted, window, preregistered.null_distribution(arrays, shifts, seed), nights)
    fitted, window = a_even(arrays)
    even = summary(fitted, window, null_even(arrays, shifts, seed), nights)
    return {"datasets": names, "measurements": sum(len(values) for values in rows.values()), "first_hjd": float(min(values[:, 0].min() for values in rows.values()) + JD_OFFSET), "last_hjd": float(max(values[:, 0].max() for values in rows.values()) + JD_OFFSET), "nights": nights, "odd": odd, "even": even, "shifts": shifts, "seed": seed, "before": before}
