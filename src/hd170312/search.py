"""The box search over the nine-sector stack, with transitgpu on the CPU.

The stack is the PDCSAP flux of every sector (cadences of QUALITY 0), divided by its biweight trend (window 0.75 d, `transitgpu.detrend`), over 0.5 to 14 d. The strongest peak is the one that the search of the
Note found; its SDE is the z-score of the box statistic over the trial periods. The Note's per-star threshold comes from the discovery pipeline's own null and is not computed here."""
import numpy as np
from transitgpu import engine
from transitgpu.detrend import biweight

from hd170312 import data

WINDOW_DAYS = 0.75
PERIOD_MIN, PERIOD_MAX = 0.5, 14.0


def stack(loader=data.load_sector, sectors=data.SECTORS):
    """(time, flux): the detrended light curves of the sectors, concatenated."""
    times, fluxes = [], []
    for sector in sectors:
        columns = loader(sector)
        times.append(columns["time"])
        fluxes.append(columns["pdcsap"] / biweight(columns["time"], columns["pdcsap"], WINDOW_DAYS, device="cpu"))
    time, flux = np.concatenate(times), np.concatenate(fluxes)
    finite = np.isfinite(flux)
    return time[finite], flux[finite]


def run(time, flux, period_min=PERIOD_MIN, period_max=PERIOD_MAX):
    """The strongest peak of the box search and the spacing of the period grid there."""
    found = engine.search(time, flux, period_min, period_max, spectrum=True, device="cpu")
    periods = np.asarray(found["periods"])
    nearest = int(np.argmin(np.abs(periods - found["period_days"])))
    step = float(np.abs(np.diff(periods)[max(nearest - 1, 0) : nearest + 1]).mean())
    return {"period_days": found["period_days"], "sde": found["sde"], "depth_ppm": found["depth_ppm"], "duration_hours": found["duration_hr"], "grid_step_days": step,
            "points": int(time.size), "period_range_searched": [found["period_min_searched_days"], found["period_max_searched_days"]]}
