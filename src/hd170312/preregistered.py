"""The test of 68.95 d against 137.90 d on the ASAS-3 nights that only 68.95 d predicts.

The model is a cos^2 dip of the TESS duration with the depth in magnitudes for each kind, times an amplitude A, at T_kind + k P with P = DeltaT / 26; epochs of even k are those of 137.90 d. The nights are those within
half a duration of an odd-k epoch and of no even-k epoch. A weighted least squares gives A with one offset per dataset. The null shifts the ephemeris by a random offset (20,000 shifts, seed 1, excluding shifts
within 1.2 d of the true phase); p = (1 + shifts with an amplitude at least as large) / (1 + 20,000). The rule fixed before the fit: 68.95 d is favoured if p < 0.01 and 0.5 < A < 1.5."""
import numpy as np

from hd170312.asas import nightly_means, read_rows, stacked

P26 = 1792.694 / 26
KINDS = [(2377.1543, 0.0855, 0.747, "shallow"), (4131.0663, 0.1536, 0.653, "deep")]  # (centre BTJD, flux depth, duration d, kind) of the Note
SHIFTS = 20_000
SEED = 1
EXCLUSION = 1.2  # d, from the true phase of either kind


def template(t, centre, duration, depth):
    """The magnitude dip (fainter is positive): -2.5 log10(1 - depth) cos^2(pi x / duration) for |x| < duration / 2."""
    x = t - centre
    out = np.zeros_like(t)
    inside = abs(x) < duration / 2
    out[inside] = -2.5 * np.log10(1 - depth) * np.cos(np.pi * x[inside] / duration) ** 2
    return out


def build(t, shift, selection, dt=0.0):
    """The template column at the epochs of 68.95 d shifted by `shift`: selection 'odd' or 'even' cycles of the ephemeris, or 'all'."""
    column = np.zeros_like(t)
    for centre, depth, duration, _ in KINDS:
        k = np.round((t - centre - shift) / P26)
        values = template(t, centre + shift + k * P26 + dt, duration, depth)
        if selection == "odd":
            values = np.where(k % 2 != 0, values, 0)
        elif selection == "even":
            values = np.where(k % 2 == 0, values, 0)
        column += values
    return column


def fit(t, m, e, datasets, column):
    """(A, its formal error, chi-square, offsets) of the weighted least squares of the magnitudes on one offset for each dataset and `column`."""
    sets = sorted(set(datasets))
    design = np.column_stack([(datasets == d).astype(float) for d in sets] + [column])
    weight = 1 / e
    a, b = design * weight[:, None], m * weight
    beta, *_ = np.linalg.lstsq(a, b, rcond=None)
    covariance = np.linalg.pinv(a.T @ a)
    return beta[-1], np.sqrt(covariance[-1, -1]), ((b - a @ beta) ** 2).sum(), beta[:-1]


def only_window(t, shift=0.0):
    """The nights within half a duration of an odd-cycle epoch and not within half a duration of an even-cycle epoch."""
    odd, even = np.zeros(len(t), bool), np.zeros(len(t), bool)
    for centre, _, duration, _ in KINDS:
        k = np.round((t - centre - shift) / P26)
        x = abs(t - (centre + shift + k * P26))
        odd |= (x < duration / 2) & (k % 2 != 0)
        even |= (x < duration / 2) & (k % 2 == 0)
    return odd & ~even


def a_only(arrays, shift=0.0):
    """(A, error, chi-square, offsets), the window of the nights tested."""
    t, m, e, datasets, _ = arrays
    window = only_window(t, shift)
    column = np.where(window, build(t, shift, "odd"), 0.0)
    return fit(t, m, e, datasets, column), window


def null_distribution(arrays, shifts=SHIFTS, seed=SEED, exclusion=EXCLUSION):
    """The amplitudes of `shifts` shifted ephemerides, drawn from a generator seeded with `seed`, excluding shifts within `exclusion` d (EXCLUSION) of the true phase of either kind; with 0.0 no shift is excluded."""
    generator = np.random.default_rng(seed)
    separation = (KINDS[1][0] - KINDS[0][0]) % P26
    values = []
    centres = (0.0, separation, (-separation) % P26)
    while len(values) < shifts:
        u = generator.uniform(0, P26)
        if any(min((u - c) % P26, (c - u) % P26) < exclusion for c in centres):
            continue
        values.append(a_only(arrays, u)[0][0])
    return np.array(values)


def decision(p, amplitude):
    """The rule fixed in advance."""
    if p < 0.01 and 0.5 <= amplitude <= 1.5:
        return "favours 68.95 d"
    return "no preference" if p > 0.10 else "weak preference"


def run(path=None, shifts=SHIFTS, seed=SEED):
    """The result of the test on the archive file: the nights, A and its error, the null, p and the decision."""
    arrays = stacked(nightly_means(read_rows(path)))
    (amplitude, sigma, _, _), window = a_only(arrays)
    null = null_distribution(arrays, shifts, seed)
    p = (1 + int((null >= amplitude).sum())) / (1 + len(null))
    unrestricted = null_distribution(arrays, shifts, seed, exclusion=0.0)  # the same test with every shift kept, including those that land on the true eclipses
    t = arrays[0]
    shallow = sum(1 for i in np.flatnonzero(window) if abs(t[i] - (KINDS[0][0] + np.round((t[i] - KINDS[0][0]) / P26) * P26)) < KINDS[0][2] / 2)
    return {"nights": int(window.sum()), "shallow_nights": shallow, "deep_nights": int(window.sum()) - shallow, "amplitude": float(amplitude), "formal_error": float(sigma), "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "template_depth_mag": [float(-2.5 * np.log10(1 - depth)) for _, depth, _, _ in KINDS], "shifts": len(null), "null_amplitudes": [float(value) for value in null], "seed": seed, "shifts_at_least_as_large": int((null >= amplitude).sum()), "p": float(p), "p_no_exclusion": float((1 + int((unrestricted >= amplitude).sum())) / (1 + len(unrestricted))), "decision": decision(p, amplitude)}
