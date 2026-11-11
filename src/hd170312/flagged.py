"""The eclipse of sector 93 and the egress of sector 101 in the SAP flux, where the cadences are flagged for scattered light (QUALITY bit 4096) and PDCSAP is empty.

The eclipse is the deep one: the SAP template of sector 102 (its depth, duration and ingress fraction) is fitted to the SAP flux of sector 93 from BTJD 3854.5 with the depth and the mid-time free. The egress is the shallow eclipse at the time that only 68.95 d
predicts: the shallow template (the mean of sectors 39 and 103) with the time fixed and the depth free, on the SAP flux of sector 101 from its first finite cadence, which lies 0.2 d after the predicted mid-time, so that only the egress is in the data."""
import numpy as np
from scipy.optimize import least_squares

from hd170312.eclipses import _model, fit_event_time

BASELINE = (3853.0, 3854.5)  # BTJD of sector 93 before the eclipse, the stretch whose median is the flux level of the fit
WINDOW_S93 = 3854.5  # the fit of the sector 93 eclipse starts here
WINDOW_S101 = (4101.11, 4102.4)  # the fit of the egress, from the first finite SAP cadence of sector 101
LEVEL_S101 = 4101.6  # the flux level of the egress is the median from here: 0.7 d after the predicted mid-time, outside the eclipse
PERIODS_26 = 26  # 68.95 d is the separation of the two shallow eclipses over 26
GUESS = {"deep": 4131.07, "shallow": (2377.16, 4169.86)}
MINIMUM_CADENCES = 100


def fit_template(time, flux, shape, centre, start, fixed_time=None):
    """The trapezoid of a fixed shape (duration, ingress fraction) over a quadratic baseline fitted to `flux` (relative): {time, depth, chisq_over_dof, flat_gain}; the time is free from `start` = (time, depth) unless `fixed_time` is given."""
    time, flux = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    if time.size < MINIMUM_CADENCES:
        raise ValueError(f"{time.size} cadences, fewer than the {MINIMUM_CADENCES} that a fit of an eclipse needs")

    def model(p):
        t0, depth = (fixed_time, p[0]) if fixed_time is not None else (p[0], p[1])
        coefficients = p[1:] if fixed_time is not None else p[2:]
        return _model([t0, depth, *shape, *coefficients], time, centre, 1.0)

    result = least_squares(lambda p: model(p) - flux, [start[1], 1.0, 0.0, 0.0] if fixed_time is not None else [start[0], start[1], 1.0, 0.0, 0.0], x_scale="jac")
    flat = np.polyfit(time - centre, flux, 2)
    far = np.abs(time - (fixed_time if fixed_time is not None else result.x[0])) >= 0.7
    sigma = float(np.std(flux[far] - np.polyval(np.polyfit(time[far] - centre, flux[far], 2), time[far] - centre)))
    chisq = float(np.sum((result.fun / sigma) ** 2))
    flat_chisq = float(np.sum(((np.polyval(flat, time - centre) - flux) / sigma) ** 2))
    parameters = len(result.x)
    return {"time": float(fixed_time if fixed_time is not None else result.x[0]), "depth": float(result.x[0] if fixed_time is not None else result.x[1]), "chisq_over_dof": chisq / (time.size - parameters), "gain_over_flat": float(1 - chisq / flat_chisq)}


def template_of(sap, guess):
    """The shape and depth of the eclipse near `guess` in a sector's SAP flux of QUALITY 0: {depth, duration, ingress_fraction}."""
    good = sap["quality"] == 0
    fit = fit_event_time(sap["time"][good], sap["sap"][good], guess)
    return {key: fit[key] for key in ("depth", "duration", "ingress_fraction")}


def fourth_eclipse(loader):
    """The deep eclipse in the flagged cadences of sector 93: {time, depth, depth_ratio (to the SAP depth of sector 102), flagged_cadences, template}; `loader(sector)` gives {time, sap, quality}."""
    template = template_of(loader(102), GUESS["deep"])
    sector = loader(93)
    stretch = (sector["time"] > BASELINE[0]) & (sector["time"] < BASELINE[1]) & (sector["quality"] == 0)
    window = sector["time"] >= WINDOW_S93
    flux = sector["sap"][window] / np.median(sector["sap"][stretch])
    found = fit_template(sector["time"][window], flux, (template["duration"], template["ingress_fraction"]), centre=3855.0, start=(3855.25, 0.15))
    return {"time": found["time"], "depth": found["depth"], "depth_ratio": found["depth"] / template["depth"], "flagged_cadences": int(np.sum(sector["quality"][window] == 4096)), "chisq_over_dof": found["chisq_over_dof"], "template": template}


def egress(loader, events):
    """The egress of the shallow eclipse at the epoch that only 68.95 d predicts, in sector 101: {predicted_time, depth, depth_ratio (to the mean shallow depth of sectors 39 and 103), gain_over_flat, first_cadence, template}; `events` are the fitted complete eclipses (hd170312.events)."""
    shallow = [template_of(loader(sector), guess) for sector, guess in ((39, GUESS["shallow"][0]), (103, GUESS["shallow"][1]))]
    template = {key: float(np.mean([t[key] for t in shallow])) for key in shallow[0]}
    predicted = events["S103"]["time"] - (events["S103"]["time"] - events["S39"]["time"]) / PERIODS_26
    sector = loader(101)
    window = (sector["time"] >= WINDOW_S101[0]) & (sector["time"] <= WINDOW_S101[1])
    time = sector["time"][window]
    flux = sector["sap"][window] / np.median(sector["sap"][window][time >= LEVEL_S101])
    found = fit_template(time, flux, (template["duration"], template["ingress_fraction"]), centre=4101.0, start=(predicted, 0.05), fixed_time=predicted)
    return {"predicted_time": float(predicted), "depth": found["depth"], "depth_ratio": found["depth"] / template["depth"], "gain_over_flat": found["gain_over_flat"], "chisq_over_dof": found["chisq_over_dof"], "first_cadence": float(time.min()), "template": template}
