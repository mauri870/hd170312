"""The mid-time of an eclipse with its error, and a weighted linear ephemeris of a few eclipses.

`fit_event_time` fits a symmetric trapezoid over a quadratic baseline. Its error has a statistical part (the covariance of the fit with the residuals' variance, times the red-noise factor beta from the scatter of
binned residuals, Winn et al. 2008, ApJ 683, 1076) and a systematic part (the scatter of the time over other baselines and windows), each times a scale that `hd170312 calibrate` measures. `fit_ephemeris` is the
weighted least squares of time = epoch(kind) + period * cycle with one epoch per kind and one period."""

import numpy as np
from scipy.optimize import least_squares

# The two scales, measured on injected eclipses (seeds 0 to 299): the statistical one is the root mean square of (time - truth) over the statistical error on dips of
# the assumed shape with white noise; the systematic one puts the share of the injected population within one total error nearest to 0.6827, the share of a normal distribution within one standard deviation.
STATISTICAL_SCALE = 0.966
SYSTEMATIC_SCALE = 1.35
BETA_BINS = (5, 10, 20)  # cadences of the bins of the time-averaged residuals
MINIMUM_BINS = 8
NOISE_LEVELS = 5.0  # the depth of the dip over its error, under which there is no dip


def _model(parameters, time, centre, scale):
    """The flux: a polynomial baseline in u = (time - centre) / scale times one minus a trapezoid of (t0, depth, duration, fraction), the ingress being `fraction` of the half duration."""
    t0, depth, duration, fraction = parameters[:4]
    half = duration / 2
    dip = depth * np.clip((half - np.abs(time - t0)) / (fraction * half), 0.0, 1.0)
    u = (time - centre) / scale
    baseline = sum(coefficient * u**power for power, coefficient in enumerate(parameters[4:]))
    return baseline * (1.0 - dip)


def _preliminary(time, flux, guess, window):
    """(depth, centre, full width at half depth, level at the guess) of the dip, from the medians of bins of 0.05 d divided by a quadratic fitted to the bins in the outer 40 per cent of the window; ValueError when
    there is none or too little data."""
    edges = np.arange(guess - window, guess + window + 0.05, 0.05)
    index = np.digitize(time, edges)
    centres, levels = [], []
    for number in range(1, edges.size):
        selected = flux[index == number]
        if selected.size >= 3:
            centres.append(0.5 * (edges[number - 1] + edges[number]))
            levels.append(np.median(selected))
    centres, levels = np.array(centres), np.array(levels)
    if centres.size < MINIMUM_BINS:
        raise ValueError(f"too few data for an event: {centres.size} bins of 0.05 d around {guess}")
    outer = np.abs(centres - guess) > 0.6 * window
    if outer.sum() < 6:
        raise ValueError(f"too few data in the outer parts of the window around {guess}: {outer.sum()} bins")
    coefficients = np.polyfit((centres[outer] - guess) / window, levels[outer], 2)
    baseline = np.polyval(coefficients, (centres - guess) / window)
    ratio = levels / baseline
    near = np.abs(centres - guess) <= 0.6 * window
    low = int(np.argmin(np.where(near, ratio, np.inf)))
    depth = 1.0 - ratio[low]
    noise = np.std(np.diff(ratio[outer])) / np.sqrt(2) + 1e-12
    if depth < NOISE_LEVELS * noise:
        raise ValueError(f"no dip near {guess}: the deepest 0.05 d bin is {depth:.2e} below its baseline, within {NOISE_LEVELS} times the noise {noise:.1e}")
    half = 1.0 - depth / 2
    left, right = low, low
    while left > 0 and ratio[left] < half:
        left -= 1
    while right < ratio.size - 1 and ratio[right] < half:
        right += 1
    return float(depth), float(centres[low]), float(centres[right] - centres[left]) + 0.05, float(np.polyval(coefficients, 0.0))


def _red_noise_factor(residual):
    """beta: the median over BETA_BINS of the scatter of the means of consecutive bins of the residuals over the white-noise expectation, at least 1 (Winn et al. 2008, ApJ 683, 1076, the time-averaging method)."""
    sigma = np.std(residual)
    ratios = []
    for size in BETA_BINS:
        count = residual.size // size
        if count < MINIMUM_BINS:
            continue
        means = residual[: count * size].reshape(count, size).mean(axis=1)
        expected = sigma / np.sqrt(size) * np.sqrt(count / (count - 1))
        ratios.append(np.std(means) / expected)
    return max(1.0, float(np.median(ratios))) if ratios else 1.0


def _fit(time, flux, guess, window, order, preliminary):
    depth0, centre0, width0, level0 = preliminary
    best = None
    for fraction in (0.3, 1.0):
        start = [centre0, depth0, 1.3 * width0, fraction] + [level0] + [0.0] * order
        lower = [centre0 - width0, 0.0, 0.3 * width0, 0.02] + [-np.inf] * (order + 1)
        upper = [centre0 + width0, 2.0 * depth0, 4.0 * width0, 1.0] + [np.inf] * (order + 1)
        x_scale = [width0 / 20, depth0 / 10, width0 / 10, 0.2] + [0.01] * (order + 1)
        result = least_squares(lambda p: _model(p, time, guess, window) - flux, start, bounds=(lower, upper), x_scale=x_scale)
        if best is None or result.cost < best.cost:
            best = result
    return best


def _complete(time, centre, width):
    """Data on both sides of the dip, beyond one full width at half depth from its centre: an ingress cut by the end of the data has its centre unobserved."""
    return time.min() < centre - width and time.max() > centre + width


def fit_event_time(time, flux, guess, half_window=1.5):
    """The mid-time of the event near `guess` in the flux within `half_window` days of it, with its error, depth (share of the baseline) and duration.

    Raises ValueError when the window holds no dip, too few data, or a partial event."""
    time, flux = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    keep = np.isfinite(time) & np.isfinite(flux) & (np.abs(time - guess) <= half_window)
    time, flux = time[keep], flux[keep]
    if time.size < 100:
        raise ValueError(f"too few data for an event: {time.size} points within {half_window} d of {guess}")
    flux = flux / np.median(flux)
    preliminary = _preliminary(time, flux, guess, half_window)
    if not _complete(time, preliminary[1], preliminary[2]):
        raise ValueError(f"the event near {guess} is partial: the data end within one width of its centre {preliminary[1]:.3f}, so its mid-time is not observed")
    main = _fit(time, flux, guess, half_window, 2, preliminary)
    parameters = main.x
    width = parameters[2] * (1.0 - parameters[3] / 2.0)  # the full width at half depth of the fitted trapezoid
    outside = np.abs(time - parameters[0]) > width
    if outside.sum() < 100:
        raise ValueError(f"the window of {half_window} d around {guess} holds {outside.sum()} points outside one width ({width:.2f} d) of the dip; use a longer window")
    noise = main.fun[outside]
    beta = _red_noise_factor(noise)
    covariance = np.linalg.pinv(main.jac.T @ main.jac)
    statistical = STATISTICAL_SCALE * float(np.std(noise) * np.sqrt(covariance[0, 0])) * beta
    estimates = [parameters[0]]
    for order, window in ((3, half_window), (2, 2 * half_window / 3), (3, 2 * half_window / 3)):
        subset = np.abs(time - guess) <= window
        estimates.append(_fit(time[subset], flux[subset], guess, window, order, preliminary).x[0])
    systematic = SYSTEMATIC_SCALE * float(np.std(estimates, ddof=1))
    return {
        "time": float(parameters[0]),
        "error": float(np.hypot(statistical, systematic)),
        "statistical_error": statistical,
        "systematic_error": systematic,
        "beta": beta,
        "depth": float(parameters[1]),
        "duration": float(parameters[2]),
        "ingress_fraction": float(parameters[3]),
        "points": int(time.size),
        "rms": float(np.std(noise)),
    }


def _checked(events):
    errors = np.array([event["error"] for event in events], dtype=float)
    if not np.all(np.isfinite(errors)) or np.any(errors <= 0):
        raise ValueError(f"every event needs a finite positive error; got {errors.tolist()}")
    kinds = list(dict.fromkeys(event["kind"] for event in events))
    if len(events) < 1 + len(kinds):
        raise ValueError(f"{len(events)} events cannot fix {1 + len(kinds)} parameters (one period and one epoch for each of {len(kinds)} kinds)")
    return errors, kinds


def fit_ephemeris(events):
    """The weighted least-squares ephemeris of events [{time, error, cycle, kind}]: time = epochs[kind] + period * cycle, with errors that the chi-square does not scale."""
    errors, kinds = _checked(events)
    times = np.array([event["time"] for event in events], dtype=float)
    cycles = np.array([event["cycle"] for event in events], dtype=float)
    design = np.zeros((len(events), 1 + len(kinds)))
    design[:, 0] = cycles
    for row, event in enumerate(events):
        design[row, 1 + kinds.index(event["kind"])] = 1.0
    weighted = design / errors[:, None]
    if np.linalg.matrix_rank(weighted) < weighted.shape[1]:
        raise ValueError("the period is not fixed: no kind has two events at different cycles")
    solution = np.linalg.lstsq(weighted, times / errors, rcond=None)[0]
    covariance = np.linalg.inv(weighted.T @ weighted)
    chi2 = float(np.sum(((times - design @ solution) / errors) ** 2))
    dof = len(events) - design.shape[1]
    return {
        "period": float(solution[0]),
        "period_error": float(np.sqrt(covariance[0, 0])),
        "epochs": {kind: float(solution[1 + index]) for index, kind in enumerate(kinds)},
        "epoch_errors": {kind: float(np.sqrt(covariance[1 + index, 1 + index])) for index, kind in enumerate(kinds)},
        "kinds": kinds,
        "covariance": covariance,
        "chi2": chi2,
        "dof": dof,
        "reduced_chi2": chi2 / dof if dof > 0 else None,
    }


def predict(fit, kind, cycle):
    """The time of the event of `kind` at `cycle` and its error, propagated through the covariance of the period and the epoch of the kind."""
    vector = np.zeros(1 + len(fit["kinds"]))
    vector[0] = cycle
    vector[1 + fit["kinds"].index(kind)] = 1.0
    return float(fit["epochs"][kind] + cycle * fit["period"]), float(np.sqrt(vector @ fit["covariance"] @ vector))


def assign_cycles(times, kinds, period):
    """The cycle number of each event within its kind, counted from the first event of the kind in the list, for a trial period."""
    first = {}
    cycles = []
    for time, kind in zip(times, kinds):
        first.setdefault(kind, time)
        cycles.append(round((time - first[kind]) / period))
    return cycles
