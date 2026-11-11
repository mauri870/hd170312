"""The three complete eclipses in the SAP flux of the MIT Quick-Look Pipeline (QLP), beside the SPOC PDCSAP fits.

QLP takes the star from the full-frame images with its own aperture and does not use the SPOC PDC (https://archive.stsci.edu/hlsp/qlp). The comparison rule, fixed before the files were fitted: depths within
DEPTH_LIMIT percentage points and mid-times within SIGMA_LIMIT combined errors of the SPOC fits (the mid-time only in sector 39: QLP times of sectors 99 to 104 carry a barycentric error). The result is reported, not judged."""
import numpy as np

from hd170312 import data, figure
from hd170312.eclipses import fit_event_time
from hd170312.events import EVENTS, fit_events

DEPTH_LIMIT = 1.5  # percentage points
SIGMA_LIMIT = 3.0


def baseline(time, flux, deepest):
    """The range of the 1 h medians 11 h to 30 h either side of `deepest` and the rise of their mean across the eclipse, relative to the local level."""
    hours, relative = figure.profile_of(np.asarray(time), np.asarray(flux), deepest)
    low, high = figure.baseline_range_of(hours, relative)
    return {"low": low, "high": high, "ramp": figure.baseline_slope_of(hours, relative)}


def run(directory=data.DATA):
    spoc = fit_events(directory)
    comparison = {}
    for name, sector, guess, _ in EVENTS:
        columns = data.load_qlp(sector, directory=directory)
        fit = fit_event_time(columns["time"], columns["flux"], guess)
        combined = float(np.hypot(fit["error"], spoc[name]["error"]))
        comparison[name] = {
            "depth": fit["depth"],
            "depth_difference_pp": abs(fit["depth"] - spoc[name]["depth"]) * 100,
            "time": fit["time"],
            "time_difference_over_error": abs(fit["time"] - spoc[name]["time"]) / combined,
        }
    deep = spoc["S102"]["time"]
    pdcsap, sap = data.load_sector(102, directory), data.load_qlp(102, directory=directory)
    ramps = baseline(pdcsap["time"], pdcsap["pdcsap"], deep), baseline(sap["time"], sap["flux"], deep)
    spoc_sap = baseline(pdcsap["time"], pdcsap["sap"], deep)
    depths = ", ".join(f"{row['depth_difference_pp']:.2f}" for row in comparison.values())
    summary = (f"QLP SAP against SPOC PDCSAP in sectors 39, 102, 103: depths differ by {depths} percentage points (limit {DEPTH_LIMIT}), the mid-time of sector 39 by {comparison['S39']['time_difference_over_error']:.2f} combined errors (limit {SIGMA_LIMIT}); "
               "the mid-times of sectors 102 and 103 are not compared, since QLP times from the second half of sector 99 to the end of sector 104 carry a barycentric error (https://archive.stsci.edu/hlsp/qlp); "
               f"the sector 102 baseline rises by {ramps[1]['ramp']:.4f} in QLP and {ramps[0]['ramp']:.4f} in PDCSAP")
    return {"summary": summary, "largest_depth_difference_pp": max(row["depth_difference_pp"] for row in comparison.values()),
            "sector_39_time_difference_over_error": comparison["S39"]["time_difference_over_error"], "sector_102_ramp": {"qlp": abs(ramps[1]["ramp"]), "spoc_sap": abs(spoc_sap["ramp"]), "pdcsap": ramps[0]["ramp"]}, "comparison": comparison, "baseline": {"pdcsap": ramps[0], "qlp": ramps[1]}}
