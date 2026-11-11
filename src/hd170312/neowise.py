"""The NEOWISE-R single-exposure (L1b) photometry of HD 170312 (https://irsa.ipac.caltech.edu/ibe/docs/neowiser/p1bs_psd/; the Explanatory Supplement, https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/): W1 and W2 profile-fit magnitudes of each exposure.

An exposure is kept with qual_frame > 0 (a frameset of good quality), cc_flags 0000 and ph_qual AA (no artifact, signal-to-noise above 10 in W1 and W2), and a position within 1 arcsec of the star. A visit is a run of exposures with no gap of
a day. The star (W1 about 7.5, W2 about 7.6) is at the saturation limits of the profile-fit photometry (W1 about 8.0, W2 about 7.5: https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec2_2.html), where the flux is overestimated by an amount
that grows with brightness. `mjd` is the mid-time of the frameset in UTC; BTJD is BJD(TDB) - 2457000 with the light-travel time to the barycentre computed for the geocentre."""
import csv

import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

from hd170312 import preregistered
from hd170312.data import DATA, STAR, verified

P26 = preregistered.P26
TARGET_FILE = "neowise/l1b_10arcsec.csv"
FIELD_FILE = "neowise/l1b_visit_field_600arcsec.csv"
GAP_DAYS = 1.0  # a gap this long between two exposures separates two visits
MAXIMUM_ARCSEC = 1.0
COMPARISON_W2 = (6.0, 11.0)  # AllWISE W2 of the comparison stars
BLEND_ARCSEC = 15.0  # a source this close to the star shares its profile
DIP_EXPOSURES = 8  # the exposures of a run in the null of the visits: the 8 exposures of the eclipse visit in the window
DIP_MINIMUM_EXPOSURES = 12  # a visit with fewer exposures is not tested
MINIMUM_FRAMES = (5, 3)  # frames of a comparison star inside the window and outside it


def distance_to_epoch(t):
    """Days from the nearest epoch of 68.95 d of either kind."""
    best = np.full(len(t), np.inf)
    for centre, _, _, _ in preregistered.KINDS:
        best = np.minimum(best, np.abs(t - (centre + np.round((t - centre) / P26) * P26)))
    return best


def btjd_of(mjd):
    """BJD(TDB) - 2457000 of UTC mjd values, for the star."""
    times = Time(mjd, format="mjd", scale="utc")
    delay = times.light_travel_time(SkyCoord(*STAR, unit="deg"), kind="barycentric", location=EarthLocation.from_geocentric(0, 0, 0, unit="m"))
    return times.tdb.jd + delay.to("d").value - 2457000.0


def separation_arcsec(ra, dec):
    a, b, c, d = map(np.radians, (STAR[0], STAR[1], ra, dec))
    h = np.sin((d - b) / 2) ** 2 + np.cos(b) * np.cos(d) * np.sin((c - a) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(h))) * 3600


def read_csv(path):
    """The columns of a TAP result as {name: array}: text for cc_flags, ph_qual and allwise_cntr, floats for the others (an empty field is NaN)."""
    with open(verified(path)) as handle:
        rows = list(csv.DictReader(handle))
    text = {"cc_flags", "ph_qual", "allwise_cntr"}
    return {name: np.array([row[name] for row in rows]) if name in text else np.array([float(row[name]) if row[name] else np.nan for row in rows]) for name in rows[0]}


def read_exposures(path=None):
    """{btjd, mjd, w1, w2, visit} of the exposures that pass the cuts, in time order; `visit` numbers the visits from 0."""
    table = read_csv(path or DATA / TARGET_FILE)
    keep = (table["qual_frame"] > 0) & (table["cc_flags"] == "0000") & (table["ph_qual"] == "AA") & (separation_arcsec(table["ra"], table["dec"]) <= MAXIMUM_ARCSEC)
    mjd = table["mjd"][keep]
    order = np.argsort(mjd)
    mjd = mjd[order]
    return {"btjd": btjd_of(mjd), "mjd": mjd, "w1": table["w1mpro"][keep][order], "w2": table["w2mpro"][keep][order], "visit": np.concatenate([[0], np.cumsum(np.diff(mjd) > GAP_DAYS)])}


def window_visit(rows):
    """(visit number, mask of its exposures in the window of the eclipse that only 68.95 d predicts) for the visit that holds that window."""
    window = preregistered.only_window(rows["btjd"])
    visits = set(rows["visit"][window])
    if len(visits) != 1:
        raise ValueError(f"the window of the eclipse that only 68.95 d predicts is in {len(visits)} visits, not one")
    return visits.pop(), window


def predicted_time(rows, visit, window):
    """The epoch of 68.95 d, of an odd cycle, that the window exposures of the visit surround, BTJD."""
    t = rows["btjd"][window]
    for centre, _, duration, _ in preregistered.KINDS:
        k = np.round((t - centre) / P26)
        if np.all(np.abs(t - (centre + k * P26)) < duration / 2):
            return float(centre + k[0] * P26)
    raise ValueError("the window exposures do not surround one epoch")


def visit_numbers(rows):
    """The numbers of the visit with the window: its exposures in and out of the window, the epoch it predicts, the time of its faintest W2 exposure and that time's distance from the epoch, and the mean magnitude of the window exposures minus the mean of the others
    in each band (positive is fainter in the window)."""
    visit, window = window_visit(rows)
    mine = rows["visit"] == visit
    other = mine & ~window
    faintest = int(np.argmax(np.where(window, rows["w2"], -np.inf)))
    epoch = predicted_time(rows, visit, window)
    return {"visits": int(rows["visit"].max()) + 1, "exposures": int(rows["btjd"].size), "window_exposures": int(window.sum()), "other_exposures": int(other.sum()), "predicted_time": epoch, "w2_minimum_time": float(rows["btjd"][faintest]),
            "w2_minimum_offset_days": float(abs(rows["btjd"][faintest] - epoch)), "w2_window_minus_other": float(rows["w2"][window].mean() - rows["w2"][other].mean()), "w1_window_minus_other": float(rows["w1"][window].mean() - rows["w1"][other].mean())}


def comparison_shifts(rows, path=None):
    """{allwise_cntr: W2 shift} of the comparison stars in the frames of the visit with the window: the mean W2 of the star's exposures in the window frames minus its mean in the other frames of the visit. A comparison star has an AllWISE counterpart with W2 from 6 to 11,
    ph_qual A in W2, qual_frame > 0 and cc_flags 0000, lies more than 15 arcsec from the star, and has at least 5 exposures in the window frames and 3 in the others."""
    visit, window = window_visit(rows)
    frames = np.round(rows["mjd"][rows["visit"] == visit], 6)
    inside = np.round(rows["mjd"][window], 6)
    table = read_csv(path or DATA / FIELD_FILE)
    keep = np.isin(np.round(table["mjd"], 6), frames) & (table["qual_frame"] > 0) & (table["cc_flags"] == "0000") & (table["allwise_cntr"] != "") & (np.array([p[1:2] for p in table["ph_qual"]]) == "A")
    low, high = COMPARISON_W2
    keep &= (table["w2mpro_allwise"] >= low) & (table["w2mpro_allwise"] <= high)
    shifts = {}
    for name in np.unique(table["allwise_cntr"][keep]):
        mine = keep & (table["allwise_cntr"] == name)
        if np.median(separation_arcsec(table["ra"][mine], table["dec"][mine])) <= BLEND_ARCSEC:
            continue
        in_window = np.isin(np.round(table["mjd"][mine], 6), inside)
        if in_window.sum() >= MINIMUM_FRAMES[0] and (~in_window).sum() >= MINIMUM_FRAMES[1]:
            shifts[str(name)] = float(table["w2mpro"][mine][in_window].mean() - table["w2mpro"][mine][~in_window].mean())
    return shifts


def window_dip(w2, width=DIP_EXPOSURES):
    """The largest, over the runs of `width` consecutive exposures of one visit, of the mean W2 of the run minus the mean of the other exposures (positive is fainter in the run)."""
    w2 = np.asarray(w2, dtype=float)
    return float(max(w2[start:start + width].mean() - np.delete(w2, np.arange(start, start + width)).mean() for start in range(len(w2) - width + 1)))


def dips_of_visits(rows):
    """{visit: window_dip of its W2} for the visits with at least DIP_MINIMUM_EXPOSURES exposures."""
    return {int(visit): window_dip(rows["w2"][rows["visit"] == visit]) for visit in np.unique(rows["visit"]) if (rows["visit"] == visit).sum() >= DIP_MINIMUM_EXPOSURES}


def dip_null(rows):
    """The dip of every visit with enough exposures, as `window_dip` gives it: the number of visits, the dip of the visit with the eclipse that only 68.95 d predicts, the largest and the median dip of the others, the median over the visits of the standard deviation of W2 in a visit,
    and the median spacing in minutes of the consecutive exposures of the eclipse visit."""
    visit, _ = window_visit(rows)
    dips = dips_of_visits(rows)
    others = np.array([dip for number, dip in dips.items() if number != visit])
    spread = [float(rows["w2"][rows["visit"] == number].std()) for number in dips]
    return {"visits_tested": len(dips), "eclipse_visit_dip": dips[visit], "largest_other_dip": float(others.max()), "median_other_dip": float(np.median(others)), "median_visit_sd": float(np.median(spread)),
            "eclipse_visit_spacing_minutes": float(np.median(np.diff(rows["mjd"][rows["visit"] == visit])) * 1440)}
