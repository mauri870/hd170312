"""The alias test: which periods P = DeltaT / N fit the three complete eclipses and predict no other, and which of them predict the deep eclipse of sector 93 too.

For each N, every eclipse the period predicts in the nine sectors is looked at in a box of 13.45 h. It is tested when half the box has cadences. A dip is present when the mean depth in the box is at least half the
observed depth of its kind and 5 standard deviations of the sector's box depths. The observed profile is injected at that time to show that a dip would have been seen. N is rejected when a tested eclipse shows no dip
although the injection is recovered."""
import itertools

import numpy as np

BOX_HOURS = 13.45
BOX = BOX_HOURS / 24  # d
CADENCES_IN_BOX = BOX_HOURS * 30  # at 30 cadences an hour
T_SHALLOW_1, T_DEEP, T_SHALLOW_2 = 2377.1543, 4131.0663, 4169.8482  # the fitted mid-times of the Note (BTJD)
SIGMA_ANCHOR = {"deep": 0.0001, "shallow": 0.0005}  # the errors of the anchors (d)
SIGMA_SEPARATION = float(np.hypot(0.0002, 0.0005))  # the error of DeltaT
EVENT_ERROR = 0.001  # d, the assumed error of the mid-time of the eclipse of sector 93: several times its fit error (0.0001 d), whose chi-square per degree of freedom is 15
EVENT_TIMES = [T_SHALLOW_1, 3855.2675, T_DEEP, T_SHALLOW_2]  # the four eclipses: the three complete ones and the deep one of sector 93, in the flagged cadences
DEEPEST = {"deep": [(4131.064, 102)], "shallow": [(2377.162, 39), (4169.843, 103)]}  # the times of the deepest points and their sectors
N_MAX = 3200  # P down to 0.56 d, the width of an eclipse
MINIMUM_SHARE = 0.5  # a tested eclipse has at least this share of its box covered
PRESENT = 0.5  # share of the observed depth of the kind that a present dip reaches
SIGMAS = 5.0


def window(time, center, half=3.0):
    """The slice of a sorted time array within `half` days of `center`."""
    return slice(np.searchsorted(time, center - half), np.searchsorted(time, center + half))


def depth_in_box(time, flux, center):
    """(depth, cadences): the deficit of the mean flux in the box against the median of 1.5 d to 3 d from the center, and the cadences in the box; depth is None with fewer than 20 cadences."""
    part = window(time, center)
    t, f = time[part], flux[part]
    distance = np.abs(t - center)
    inside = distance < BOX / 2
    if inside.sum() < 20:
        return None, int(inside.sum())
    outside = (distance > 1.5) & (distance < 3.0)
    level = np.median(f[outside]) if outside.sum() > 20 else np.median(flux)
    return float(-(f[inside].mean() / level - 1)), int(inside.sum())


def profile_of(sectors, deepest):
    """(hours, deficit): the means over the events of the deficit against hours from the deepest point, in 0.5 h medians from -30 h to +30 h."""
    edges = np.arange(-30, 30.01, 0.5)
    stacked = []
    for center, sector in deepest:
        time, flux = sectors[sector]
        part = window(time, center, 3.0)
        t, f = time[part], flux[part]
        level = np.median(f[np.abs(t - center) > 1.5])
        near = np.abs(t - center) < 30 / 24
        hours, deficit = (t[near] - center) * 24, -(f[near] / level - 1)
        stacked.append([np.median(deficit[(hours >= lo) & (hours < hi)]) if np.any((hours >= lo) & (hours < hi)) else 0.0 for lo, hi in itertools.pairwise(edges)])
    return 0.5 * (edges[1:] + edges[:-1]), np.mean(stacked, axis=0)


def sigma_of(time, flux, event_times):
    """The scatter of the box depths of a sector away from every event (3 d) and its own edges, from boxes at a spacing of one box; None with fewer than 11 boxes."""
    depths = []
    for center in np.arange(time.min() + 2, time.max() - 2, BOX):
        if min(abs(center - event) for event in event_times) < 3:
            continue
        depth, cadences = depth_in_box(time, flux, center)
        if depth is not None and cadences > 0.8 * CADENCES_IN_BOX:
            depths.append(depth)
    return float(np.std(depths)) if len(depths) > 10 else None


def prepare(sectors, deepest=DEEPEST, event_times=EVENT_TIMES):
    """The observed profile and depth of each kind and the box scatter of each sector."""
    grid = None
    profile = {}
    for kind, events in deepest.items():
        grid, profile[kind] = profile_of(sectors, events)
    observed = {kind: float(values[np.abs(grid) < BOX_HOURS / 2].mean()) for kind, values in profile.items()}
    sigma = {sector: sigma_of(time, flux, event_times) for sector, (time, flux) in sectors.items()}
    return {"grid": grid, "profile": profile, "observed_depth": observed, "sigma_box": sigma}


def test_eclipse(sectors, prepared, kind, center, sector):
    """The record of one predicted eclipse: its cadences and box share and, when tested, its depth, the depth needed, whether a dip is present and whether the injection is recovered."""
    time, flux = sectors[sector]
    depth, points = depth_in_box(time, flux, center)
    share = points / CADENCES_IN_BOX
    record = {"kind": kind, "sector": sector, "time": float(center), "cadences_in_box": points, "share_of_box": share, "tested": False}
    if depth is None or share < MINIMUM_SHARE:
        return record
    sigma = prepared["sigma_box"][sector]
    observed = prepared["observed_depth"][kind]
    needed = max(PRESENT * observed, SIGMAS * sigma) if sigma else PRESENT * observed
    part = window(time, center, 3.0)
    injected = flux[part] * (1 - np.interp((time[part] - center) * 24, prepared["grid"], prepared["profile"][kind], left=0.0, right=0.0))
    recovered, _ = depth_in_box(time[part], injected, center)
    record.update(tested=True, depth=depth, needed=needed, present=bool(depth >= needed), recovered=recovered, recovery_detected=bool(recovered is not None and (recovered - depth) >= needed))
    return record


def alias_table(sectors, prepared=None, anchors=((("deep", T_DEEP)), ("shallow", T_SHALLOW_2)), separation=T_SHALLOW_2 - T_SHALLOW_1, n_values=range(1, N_MAX + 1), event_times=EVENT_TIMES):
    """{N: summary} for each N: the period, the predicted eclipses in the data, how many were tested and how many are of the four known events, the eclipses that reject it, and whether it is rejected."""
    prepared = prepare(sectors) if prepared is None else prepared
    summary = {}
    for n_value in n_values:
        period = separation / n_value
        counts = {"eclipses_in_the_data": 0, "tested": 0, "of_which_events": 0}
        refuting = []
        for kind, anchor in anchors:
            for sector, (time, flux) in sectors.items():
                low = int(np.ceil((time.min() - BOX / 2 - anchor) / period))
                high = int(np.floor((time.max() + BOX / 2 - anchor) / period))
                for n in range(low, high + 1):
                    center = anchor + n * period
                    record = test_eclipse(sectors, prepared, kind, center, sector)
                    counts["eclipses_in_the_data"] += 1
                    counts["of_which_events"] += int(min(abs(center - event) for event in event_times) < 0.3)
                    counts["tested"] += int(record["tested"])
                    if record["tested"] and not record["present"] and record["recovery_detected"]:
                        refuting.append({"kind": kind, "n": n, "sector": sector, "time": round(float(center), 3), "share_of_box": round(record["share_of_box"], 2)})
        summary[n_value] = {"period": float(period), **counts, "rejected": bool(refuting), "refuting_eclipses": refuting}
    return summary


def fourth_event_table(n_values, anchors, separation, kind, time, error):
    """{N: {cycles, residual_days, tolerance_days, fits}}: the eclipse of `kind` that the period separation / N predicts nearest to the observed eclipse at `time`, and whether it lies within three times the combined error of the eclipse, the anchor of the kind and
    the cycles times the error of the separation over N. A period that predicts no eclipse of that kind at the observed one is rejected by it."""
    anchor = dict(anchors)[kind]
    table = {}
    for n_value in n_values:
        period = separation / n_value
        cycles = round((time - anchor) / period)
        residual = anchor + cycles * period - time
        tolerance = 3 * float(np.sqrt(error**2 + SIGMA_ANCHOR[kind] ** 2 + (cycles * SIGMA_SEPARATION / n_value) ** 2))
        table[n_value] = {"cycles": cycles, "residual_days": float(residual), "tolerance_days": tolerance, "fits": bool(abs(residual) <= tolerance)}
    return table


def survivors(summary):
    """The values of N that no eclipse rejects."""
    return [n for n, row in summary.items() if not row["rejected"]]
