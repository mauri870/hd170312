"""Figure 1 of the Note: (a) the coverage of the nine sectors and the eclipses each period predicts; (b) the four eclipses aligned on their deepest point; (c) the ASAS-3 nights that only 68.95 d predicts."""
import matplotlib

matplotlib.use("Agg")
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hd170312 import asassn, data, neowise, preregistered
from hd170312.asas import nightly_means, read_rows, stacked
from hd170312.forecast import SECTOR_114

COLOURS = {"deep": "#b2182b", "shallow": "#2166ac"}
ORANGE = "#e08214"  # the eclipses that only 68.95 d predicts: their lines, the egress of sector 101 and the eclipse of 2027 March 11
PHASE_STEP = 0.002  # the gap in phase that ends a bar of panel (a)
SECTORS = list(data.SECTORS)
# The four eclipses drawn in the SAP flux, so that the PDC tilt beside the sector 102 eclipse does not show: the time of the deepest point of each complete eclipse and the fitted mid-time of the deep eclipse of sector 93 in the SAP flux of its flagged cadences (BTJD), and their sectors.
EVENTS = [
    {"sector": 39, "time": 2377.162, "kind": "shallow", "source": "sap"},
    {"sector": 93, "time": 3855.2675, "kind": "deep", "source": "sap"},
    {"sector": 102, "time": 4131.064, "kind": "deep", "source": "sap"},
    {"sector": 103, "time": 4169.843, "kind": "shallow", "source": "sap"},
]
LABELS = {39: "S39 shallow", 93: "S93 deep", 102: "S102 deep", 103: "S103 shallow"}
DELTA_T = 1792.694  # the separation of the two shallow eclipses, the period of the text: P = DELTA_T / N
DEEP_ANCHOR, SHALLOW_ANCHOR = 4131.0663, 4169.8482  # the fitted times of the deeper and the second shallower eclipse
P26 = DELTA_T / 26
P13 = DELTA_T / 13
EGRESS = {"sector": 101, "time": SHALLOW_ANCHOR - P26, "kind": "shallow"}  # the epoch that only 68.95 d predicts, whose egress is the start of sector 101
LEVEL_S101 = 4101.6  # the flux level of the egress is the median of the SAP flux from here
_SAP = {}


def sap_data(sector):
    """{time, sap, quality} of a sector with its flagged cadences, read once (data.load_sector_sap)."""
    if sector not in _SAP:
        _SAP[sector] = data.load_sector_sap(sector)
    return _SAP[sector]


def epochs(period, anchor, low=1650.0, high=4600.0):
    first, last = int(np.ceil((low - anchor) / period)), int(np.floor((high - anchor) / period))
    return [anchor + n * period for n in range(first, last + 1)]


def profile_of(time, flux, event_time, window=1.5):
    """Hours from `event_time` and the relative flux (1.00 is the median flux 1.5 to 3 d from it, the local level) within `window` days."""
    near = np.abs(time - event_time) < 3.0
    outside = near & (np.abs(time - event_time) > window)
    level = np.median(flux[outside]) if outside.sum() > 20 else np.median(flux)
    inside = np.abs(time - event_time) < window
    return (time[inside] - event_time) * 24, flux[inside] / level


def profile(event, window=1.5):
    """(hours, relative flux, flagged) of a sector's event from its SAP flux with the cadences flagged for scattered light kept (profile_of), at the level of the stretch from 12 to 36 h before the minimum, which the eclipse of sector 93, whose data end 7 h after it, needs; `flagged` marks the cadences of QUALITY 4096."""
    sector = sap_data(event["sector"])
    keep = np.abs(sector["time"] - event["time"]) < window
    hours, relative = profile_of(sector["time"], sector["sap"], event["time"], window)
    before = (hours < -12) & (hours > -36)
    return hours, relative / np.median(relative[before]), sector["quality"][keep] == 4096


def egress_profile():
    """(hours from the epoch that only 68.95 d predicts, relative flux) of the SAP flux at the start of sector 101, whose level is the median from LEVEL_S101."""
    sector = sap_data(EGRESS["sector"])
    hours = (sector["time"] - EGRESS["time"]) * 24
    keep = hours < 36
    return hours[keep], sector["sap"][keep] / np.median(sector["sap"][sector["time"] >= LEVEL_S101])


def hourly_medians(hours, relative):
    """Centres (hours) and medians of 1 h bins from -36 h to 36 h; a bin with no point is NaN."""
    edges = np.arange(-36, 36.01, 1.0)
    centres = 0.5 * (edges[1:] + edges[:-1])
    medians = [np.median(relative[(hours >= lo) & (hours < hi)]) if np.any((hours >= lo) & (hours < hi)) else np.nan for lo, hi in itertools.pairwise(edges)]
    return centres, np.array(medians)


def baseline_range_of(hours, relative, beyond=11.0, within=30.0):
    """(lowest, highest) 1 h median of a profile between `beyond` and `within` hours from its centre, either side."""
    centres, medians = hourly_medians(hours, relative)
    side = (np.abs(centres) > beyond) & (np.abs(centres) <= within)
    return float(np.nanmin(medians[side])), float(np.nanmax(medians[side]))


def baseline_slope_of(hours, relative, beyond=11.0, within=30.0):
    """The mean of the 1 h medians between `beyond` and `within` hours after the centre minus their mean between the same hours before it: the ramp across the eclipse."""
    centres, medians = hourly_medians(hours, relative)
    after = (centres > beyond) & (centres <= within)
    before = (centres < -beyond) & (centres >= -within)
    return float(np.nanmean(medians[after]) - np.nanmean(medians[before]))


def phase(t):
    """The phase of the 137.90 d orbit (DELTA_T / 13) of a time in BTJD, 0 at the shallow eclipse, in 0 to 1."""
    return ((np.asarray(t) - SHALLOW_ANCHOR) / P13) % 1.0


DEEP_PHASE = float(phase(DEEP_ANCHOR))
# The four lines of panel (a): (phase, kind, only 68.95 d predicts it): the eclipses that both periods predict and, half an orbit later, those that only 68.95 d predicts.
LINES = [(0.0, "shallow", False), (DEEP_PHASE, "deep", False), (0.5, "shallow", True), ((DEEP_PHASE + 0.5) % 1.0, "deep", True)]
ROWS = [*SECTORS, 114]
YEARS = {13: 2019, 39: 2021, 66: 2023, 93: 2025, 94: 2025, 100: 2026, 101: 2026, 102: 2026, 103: 2026}


def row_y(sector):
    """The height of the row of a sector in panel (a), the first sector at the top."""
    return len(ROWS) - 1 - ROWS.index(sector)


def marker_phase(t):
    """The phase of a time for a marker: 0 to 1, a phase just under 1 (the eclipse of sector 103, a moment before the anchor) drawn at 0 and not at 1."""
    value = float(phase(t))
    return value - 1.0 if value > 0.99 else value


def decisive_eclipse():
    """The deep eclipse in the window of sector 114 that only 68.95 d predicts, BTJD."""
    return next(t for t in epochs(P26, DEEP_ANCHOR, low=4400, high=4600) if SECTOR_114[0] <= t <= SECTOR_114[1])


def segments(times, step=PHASE_STEP):
    """The intervals (first phase, last phase) of the cadences, one for each run of cadences with no gap longer than `step` in phase; a single cadence is an interval of zero width."""
    phases = np.sort(phase(times))
    if phases.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(phases) > step)
    starts, ends = np.concatenate([[0], breaks + 1]), np.concatenate([breaks, [phases.size - 1]])
    return [(float(phases[first]), float(phases[last])) for first, last in zip(starts, ends, strict=True)]


NEOWISE_GREY = "0.35"
SN_GREY = "0.72"  # the ASAS-SN nights, lighter than the ASAS-3 ones


def phase_diagram(ax, arrays, exposures, sn_arrays):
    """Panel (a): one row for each sector with the phases of the 137.90 d orbit that its cadences cover (dark grey QUALITY 0, light grey the cadences flagged for scattered light; sector 114 planned, hatched), the lines of the eclipses that both
    periods predict and of those that only 68.95 d predicts, the four observed eclipses, the egress of sector 101 and the eclipse of 2027 March 11, above sector 13 a row of the NEOWISE visits (`exposures`, as neowise.read_exposures gives them) at their phase, the visit
    with the eclipse that only 68.95 d predicts in orange, above it a row of the ASAS-3 nights (`arrays`, as asas.stacked gives them) at their phase, the nights of the pre-registered test in orange and taller, and above that a row of the ASAS-SN nights (`sn_arrays`) in light grey, each night in a window of the same test as a hollow orange rectangle the width of its window, centred on its epoch."""
    for sector in ROWS:
        y = row_y(sector)
        if sector == 114:
            for a, b in segments(np.arange(*SECTOR_114, 0.01)):
                ax.add_patch(plt.Rectangle((a, y - 0.32), b - a, 0.64, fc="white", ec="0.3", ls="--", lw=0.6, hatch="////", zorder=2))
            label = "TESS S114 (2027, planned)"
        else:
            columns = sap_data(sector)
            good = columns["quality"] == 0
            for mask, colour in ((good, "0.55"), (~good, "0.8")):
                for a, b in segments(columns["time"][mask]):
                    ax.add_patch(plt.Rectangle((a, y - 0.32), b - a, 0.64, color=colour, lw=0, zorder=2))
            label = f"TESS S{sector} ({YEARS[sector]})"
        ax.text(-0.012, y, label, ha="right", va="center", fontsize=6.5)
    visit, _ = neowise.window_visit(exposures)
    for number in range(int(exposures["visit"].max()) + 1):
        for a, b in segments(exposures["btjd"][exposures["visit"] == number]):
            ax.add_patch(plt.Rectangle((a, len(ROWS) - 0.32), b - a, 0.64, color=NEOWISE_GREY, lw=0, zorder=2))
    visit_phases = phase(exposures["btjd"][exposures["visit"] == visit])
    ax.plot([(visit_phases.min() + visit_phases.max()) / 2], [len(ROWS)], marker="v", ms=6, color=ORANGE, mec="k", mew=0.5, ls="", zorder=5)  # a marker, as the other eclipses seen have, so the visit is not read as an ASAS-3 tick
    ax.text(-0.012, len(ROWS), "NEOWISE (2014–24)", ha="right", va="center", fontsize=6.5)
    nights, tested = arrays[0], preregistered.only_window(arrays[0])
    row = len(ROWS) + 1  # the row of the ASAS-3 nights, above the NEOWISE visits
    ax.vlines(phase(nights[~tested]), row - 0.3, row + 0.3, color="0.55", lw=0.35, zorder=2)
    ax.vlines(phase(nights[tested]), row - 0.42, row + 0.42, color=ORANGE, lw=1.3, zorder=3)
    ax.text(-0.012, row, "ASAS-3 (2000–09)", ha="right", va="center", fontsize=6.5)
    sn_nights, sn_tested = sn_arrays[0], preregistered.only_window(sn_arrays[0])
    sn_row = row + 1
    ax.vlines(phase(sn_nights), sn_row - 0.3, sn_row + 0.3, color=SN_GREY, lw=0.35, zorder=2)
    for night in sn_nights[sn_tested]:
        centre, _, duration, _ = min(preregistered.KINDS, key=lambda kind: abs(night - (kind[0] + np.round((night - kind[0]) / preregistered.P26) * preregistered.P26)))
        epoch = centre + np.round((night - centre) / preregistered.P26) * preregistered.P26
        ax.add_patch(plt.Rectangle((phase(epoch) - duration / P13 / 2, sn_row - 0.42), duration / P13, 0.84, fill=False, ec=ORANGE, lw=0.6, zorder=3))  # the window of the test: half a duration each side of the epoch
    ax.text(-0.012, sn_row, "ASAS-SN (2014–26)", ha="right", va="center", fontsize=6.5)
    for x, kind, only in LINES:
        colour = ORANGE if only else COLOURS[kind]
        ax.plot([x, x], [-0.55, len(ROWS) + 2.6], color=colour, lw=1.4 if only else 1.0, ls=(0, (3, 1.5)) if only else "-", zorder=1)  # the line stops under its label
        ax.text(x, len(ROWS) + 2.75, ("68.95 d only\n" if only else "both periods\n") + kind, ha="center", va="bottom", fontsize=6, color=colour)
    for event in EVENTS:
        ax.plot([marker_phase(event["time"])], [row_y(event["sector"])], marker="s" if event["kind"] == "shallow" else "o", ms=6, color=COLOURS[event["kind"]], mec="k", mew=0.5, zorder=5)
    ax.plot([marker_phase(EGRESS["time"])], [row_y(EGRESS["sector"])], marker="D", ms=5.5, color=ORANGE, mec="k", mew=0.5, zorder=5)
    ax.plot([marker_phase(decisive_eclipse())], [row_y(114)], marker="o", ms=8, mfc="white", mec=ORANGE, mew=1.4, zorder=5)
    ax.set_xlim(-0.02, 1.0)
    ax.set_ylim(-0.7, len(ROWS) + 3.9)
    ax.set_yticks([])
    ax.set_xlabel("phase of the 137.90 d orbit (0 = shallow eclipse)", fontsize=7)
    ax.tick_params(labelsize=6.5)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color="0.55", label="TESS data"),
        plt.Rectangle((0, 0), 1, 1, color="0.8", label="flagged data"),
        plt.Rectangle((0, 0), 1, 1, fc="white", ec="0.3", ls="--", lw=0.6, hatch="////", label="planned"),
        plt.Line2D([], [], ls="", marker="s", ms=5, color=COLOURS["shallow"], mec="k", mew=0.5, label="shallow eclipse seen"),
        plt.Line2D([], [], ls="", marker="o", ms=5, color=COLOURS["deep"], mec="k", mew=0.5, label="deep eclipse seen"),
        plt.Line2D([], [], ls="", marker="D", ms=4.5, color=ORANGE, mec="k", mew=0.5, label="possible egress, flagged data (68.95 d only)"),
        plt.Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec=ORANGE, mew=1.2, label="predicted eclipse, 2027 March 11"),
        plt.Rectangle((0, 0), 1, 1, color=NEOWISE_GREY, label="NEOWISE visit"),
        plt.Line2D([], [], ls="", marker="v", ms=5, color=ORANGE, mec="k", mew=0.5, label="NEOWISE visit with the eclipse only 68.95 d predicts"),
        plt.Line2D([], [], ls="", marker="|", ms=6, mew=0.6, color="0.55", label="ASAS-3 night"),
        plt.Line2D([], [], ls="", marker="|", ms=9, mew=1.3, color=ORANGE, label="ASAS-3 night in an eclipse only 68.95 d predicts"),
        plt.Line2D([], [], ls="", marker="|", ms=6, mew=0.6, color=SN_GREY, label="ASAS-SN night"),
        plt.Rectangle((0, 0), 1, 1, fill=False, ec=ORANGE, lw=0.6, label="ASAS-SN night in an eclipse only 68.95 d predicts"),
    ]
    ax.legend(handles=handles, fontsize=6, frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.45, -0.2), handletextpad=0.3, columnspacing=0.9)


def eclipses(ax):
    for event in EVENTS:
        hours, relative, flagged = profile(event)
        colour = COLOURS[event["kind"]]
        ax.plot(hours[~flagged], relative[~flagged], ".", ms=1.2, color=colour, alpha=0.3, rasterized=True)
        ax.plot(hours[flagged], relative[flagged], "o", ms=1.6, mfc="none", mec=colour, mew=0.3, alpha=0.5, rasterized=True)
        centres, medians = hourly_medians(hours, relative)
        ax.plot(centres, medians, "-", lw=1.2, color=colour, label=LABELS[event["sector"]])
    hours, relative = egress_profile()
    ax.plot(hours, relative, ".", ms=1.4, color="k", alpha=0.45, rasterized=True, label="_egress")  # last, so that it stays on top of the curves; its legend entry is the handle below
    ax.axhline(1.0, color="k", lw=0.4)
    ax.set_xlim(-30, 30)
    ax.set_ylim(0.825, 1.02)
    ax.set_xlabel("hours from the deepest point", fontsize=7.5)
    ax.set_ylabel("relative flux", fontsize=7.5)
    handles, _ = ax.get_legend_handles_labels()
    handles.append(plt.Line2D([], [], ls="-", lw=1.2, color="k", label="S101, possible egress (flagged)"))  # a line like the other entries: the points are too small to show in a legend
    ax.legend(handles=handles, fontsize=6.2, frameon=False, loc="lower left", handlelength=1.2, borderaxespad=0.2, labelspacing=0.3)


def error_scale(arrays):
    """The root of the chi-square per degree of freedom of the fit of the test (one offset for each dataset and the amplitude): the factor by which the quoted errors of the nights are widened in panel (c)."""
    (_, _, chi2, _), _ = preregistered.a_only(arrays)
    return float(np.sqrt(chi2 / (len(arrays[0]) - (len(set(arrays[3])) + 1))))


def stored_null(path=None):
    """The amplitudes of the shifted ephemerides that the test wrote into out/asas-preregistered.json, or None when the file is not there or holds none: the figure draws them and does not compute them again."""
    path = Path(path or Path(__file__).resolve().parents[2] / "out" / "asas-preregistered.json")
    if not path.exists():
        return None
    values = json.loads(path.read_text()).get("null_amplitudes")
    return None if values is None else np.asarray(values, dtype=float)


def asas(ax, path=None, null=None, arrays=None):
    """Panel (c): the nights of the pre-registered test and the other nights within 2 d of the same epochs, with the two templates at unit amplitude and the error bars widened by the root of the chi-square per degree of freedom of the fit, and an inset with the histogram of the null; returns (nights, amplitude, error); the null is `null`, else the stored one, else computed."""
    arrays = stacked(nightly_means(read_rows(path))) if arrays is None else arrays
    T, M, E, DS, _ = arrays
    (amplitude, sigma, _, offsets), win = preregistered.a_only(arrays)
    scaled = E * error_scale(arrays)
    sets = sorted(set(DS))
    offset = sum((DS == d) * offsets[i] for i, d in enumerate(sets))
    print(f"nights {int(win.sum())}, A_only {amplitude:.3f} +- {sigma:.3f}")
    for centre, depth, dur, kind in preregistered.KINDS:
        k = np.round((T - centre) / P26)
        x = T - (centre + k * P26)
        odd = k % 2 != 0
        near_even = np.zeros(len(T), bool)
        for other, _, other_duration, _ in preregistered.KINDS:
            kk = np.round((T - other) / P26)
            near_even |= (abs(T - (other + kk * P26)) < other_duration / 2 + 0.2) & (kk % 2 == 0)
        context = odd & (abs(x) < 2.0) & ~win & ~near_even
        inside = odd & win & (abs(x) < dur / 2)
        ax.errorbar(x[context], (M - offset)[context], scaled[context], fmt="o", ms=2, color="0.65", ecolor="0.8", elinewidth=0.5, zorder=2)
        ax.errorbar(x[inside], (M - offset)[inside], scaled[inside], fmt="o", ms=3.5, color=COLOURS[kind], ecolor=COLOURS[kind], elinewidth=0.8, zorder=4, label=f"{kind}, {int(inside.sum())} nights")
        grid = np.linspace(-1.2, 1.2, 400)
        ax.plot(grid, preregistered.template(grid, 0.0, dur, depth), "-", color=COLOURS[kind], lw=1.2, zorder=3)
    ax.axhline(0, color="k", lw=0.4)
    ax.set_ylim(0.30, -0.10)  # fainter is down
    ax.set_xlim(-1.6, 1.6)
    ax.set_xlabel("days from the epoch that only 68.95 d predicts", fontsize=7.5)
    ax.set_ylabel("ASAS-3 V (mag)", fontsize=7.5)
    ax.legend(fontsize=7, frameon=False, loc="lower left", handlelength=1.0, borderaxespad=0.2, labelspacing=0.3)
    null = stored_null() if null is None else null
    null = preregistered.null_distribution(arrays) if null is None else null
    p = (1 + int((null >= amplitude).sum())) / (1 + null.size)
    inset = ax.inset_axes([0.60, 0.17, 0.37, 0.27])
    inset.hist(null, bins=60, color="0.6", log=True)
    inset.axvline(amplitude, color="k", lw=1)
    inset.set_yticks([])
    inset.tick_params(labelsize=5)
    inset.set_xlabel("amplitude, shifted epochs", fontsize=5, labelpad=1)
    inset.text(amplitude, inset.get_ylim()[1] * 0.3, f" {amplitude:.2f}\n p={p:.3f}", fontsize=5, va="top")
    return int(win.sum()), amplitude, sigma


def make(out, asas_path=None, null=None, neowise_path=None, asassn_path=None):
    """Draw the figure into the file `out` and into a PDF of the same name, the format the Note uses; return the numbers of panel (c)."""
    plt.rcParams.update({"font.size": 7, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7})
    fig = plt.figure(figsize=(6.5, 6.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], hspace=0.62, wspace=0.34)
    a, b, c = fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    arrays = stacked(nightly_means(read_rows(asas_path)))  # the ASAS-3 file is read once, for panels (a) and (c)
    phase_diagram(a, arrays, neowise.read_exposures(neowise_path), asassn.read_nights(asassn_path)[0])
    eclipses(b)
    numbers = asas(c, None, null, arrays)
    for axis, label in ((a, "(a)"), (b, "(b)"), (c, "(c)")):
        axis.set_title(label, fontsize=8, loc="left", fontweight="bold")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(Path(out).with_suffix(".pdf"), bbox_inches="tight")  # the AAS accepts PDF figures and not PNG; the PNG is for the README
    plt.close(fig)
    return numbers
