"""The command line: python -m hd170312 COMMAND (uv run hd170312 COMMAND)."""
import argparse
import sys

COMMANDS = {}


def command(function):
    COMMANDS[function.__name__.replace("_", "-")] = function
    return function


@command
def fetch(arguments):
    """Download the TESS light curves and the ASAS-3 file into data/ and check them against data/checksums.json (--record adds the files that have no checksum)."""
    from hd170312 import data

    paths = data.fetch(record=arguments.record)
    for name, path in paths.items():
        print(f"{name}: {data.sha256(path)} {path.stat().st_size} bytes")
    return 0


@command
def eclipses(arguments):
    """Fit the mid-time, depth and duration of the three complete eclipses (sectors 39, 102, 103) with their errors; write out/eclipses.json."""
    from hd170312.events import fit_events

    return output("eclipses", fit_events())


@command
def radius(arguments):
    """The lower bound on the radius of the body that eclipses A, for each kind of eclipse and each radius of A, with and without limb darkening; write out/radius.json."""
    from hd170312 import radius as module
    from hd170312.events import fit_events

    return output("radius", module.table(fit_events()))


@command
def flagged(arguments):
    """The deep eclipse in the flagged cadences of sector 93 (the SAP template of sector 102, depth and mid-time free) and the egress of the shallow eclipse in sector 101 (time fixed at the 68.95 d prediction); write out/flagged.json."""
    from hd170312 import data
    from hd170312 import flagged as module
    from hd170312.events import fit_events

    return output("flagged", {"fourth_eclipse": module.fourth_eclipse(data.load_sector_sap), "egress": module.egress(data.load_sector_sap, fit_events())})


@command
def alias(arguments):
    """The alias test: which periods DeltaT / N, N = 1 to 3200, fit the three complete eclipses and predict no other eclipse in the nine sectors, and which of those also predict the deep eclipse of sector 93; write out/alias.json. About a minute of one core."""
    from hd170312 import alias as module
    from hd170312 import data, flagged

    sectors = {}
    for sector in data.SECTORS:
        columns = data.load_sector(sector)
        sectors[sector] = (columns["time"], columns["pdcsap"])
    prepared = module.prepare(sectors)
    table = module.alias_table(sectors, prepared)
    fit_three = module.survivors(table)
    fourth = flagged.fourth_eclipse(data.load_sector_sap)
    anchors = (("deep", module.T_DEEP), ("shallow", module.T_SHALLOW_2))
    by_fourth = module.fourth_event_table(table, anchors, module.T_SHALLOW_2 - module.T_SHALLOW_1, "deep", fourth["time"], module.EVENT_ERROR)
    kept = [n for n in fit_three if by_fourth[n]["fits"]]
    return output("alias", {"fit_three_eclipses": fit_three, "fourth_eclipse": {"time": fourth["time"], "assumed_error": module.EVENT_ERROR, "by_N": {n: by_fourth[n] for n in range(1, 31)}}, "rule": {"minimum_share_of_box": module.MINIMUM_SHARE, "present_share_of_observed_depth": module.PRESENT, "sigmas": module.SIGMAS}, "survivors": kept, "rejected": len(table) - len(kept), "tested_values": len(table), "observed_depth": prepared["observed_depth"], "sigma_box_ppm": {s: (v * 1e6 if v else None) for s, v in prepared["sigma_box"].items()},
                           "rejected_by_one_eclipse_up_to_15": [n for n, row in table.items() if n <= 15 and row["rejected"] and len(row["refuting_eclipses"]) == 1], "by_N": {n: row for n, row in table.items() if n <= 60 or not row["rejected"]}})


@command
def forecast(arguments):
    """Fit the ephemeris of each period to the three eclipse times and predict the eclipses in TESS sector 114, with the deciding eclipse of the shortest period; write out/forecast.json."""
    from hd170312 import forecast as module
    from hd170312.events import fit_events

    found = fit_events()
    events = [{"time": e["time"], "error": e["error"], "kind": e["kind"]} for e in found.values()]
    rows = module.forecast(events)
    separation = [e for e in events if e["kind"] == "shallow"][-1]["time"] - next(e for e in events if e["kind"] == "shallow")["time"]
    return output("forecast", {"separation": separation, "separation_over_139": separation / 139, "deciding": module.deciding_eclipse(rows), "sector_114": module.SECTOR_114, "by_N": rows, "predicted_eclipses_in_the_window": module.minimum_gap(rows)})


@command
def note(arguments):
    """Build the PDF of the Note, note/main.pdf, with tectonic (https://tectonic-typesetting.github.io/; it fetches the LaTeX packages it needs on its first run)."""
    import shutil
    import subprocess
    from pathlib import Path

    directory = Path(__file__).resolve().parents[2] / "note"
    if shutil.which("tectonic") is None:
        print("tectonic is not installed", file=sys.stderr)
        return 1
    return subprocess.run(["tectonic", "main.tex"], cwd=directory, check=False).returncode


@command
def asas_preregistered(arguments):
    """The test of 68.95 d against 137.90 d on the ASAS-3 nights that only 68.95 d predicts (20,000 shifts, seed 1), in grade A in the aperture MAG_3; write out/asas-preregistered.json."""
    from hd170312 import preregistered

    return output("asas-preregistered", preregistered.run())


@command
def asassn(arguments):
    """The ASAS-SN photometry of saturated stars (the file of data/asassn, checked against its checksum): the amplitude of the deficit at the nights that only 68.95 d places in eclipse and, as a control, at the nights that both periods place in eclipse, each with 20,000 shifted ephemerides (seed 1), and the same two tests on the nights before BTJD 3900, before the drift of the cameras; write out/asassn.json."""
    from hd170312 import asassn as module

    return output("asassn", {**module.run(), "before_3900": module.run(before=3900.0)})


@command
def figure(arguments):
    """Draw Figure 1 of the Note into out/figure.png and out/figure.pdf."""
    from pathlib import Path

    from hd170312 import figure as module

    out = Path(__file__).resolve().parents[2] / "out"
    out.mkdir(exist_ok=True)
    nights, amplitude, error = module.make(out / "figure.png")
    print(f"out/figure.png, out/figure.pdf: {nights} nights, A = {amplitude:.3f} +- {error:.3f}")
    return 0


@command
def qlp(arguments):
    """Fit the three complete eclipses in the QLP SAP flux (sectors 39, 102, 103) and compare them with the SPOC fits; write out/qlp.json."""
    from hd170312 import qlp as module

    return output("qlp", module.run())


@command
def search(arguments):
    """Search the nine-sector stack (PDCSAP, biweight 0.75 d) for box transits over 0.5 to 14 d with transitgpu on the CPU; write out/search.json. About 40 s; needs `uv sync --locked --group search`."""
    from hd170312 import claims

    if not claims.search_installed():
        print(f"skipped: {claims.SEARCH_HINT}", file=sys.stderr)
        return 0
    from hd170312 import search as module

    return output("search", module.run(*module.stack()))


@command
def asas_pooled(arguments):
    """The ASAS-3 light curve (grades A and B, the aperture with the lowest median frame error) against the eclipse trains of the seven periods: the shift tests at 50,000 shifts of 68.95 d, its epochs that 137.90 d leaves out and 137.90 d; write out/asas-pooled.json. About a minute of one core."""
    from hd170312 import asas, pooled
    from hd170312.events import fit_events

    events = [{"time": e["time"], "error": e["error"], "kind": e["kind"]} for e in fit_events().values()]
    time, flux, _, _ = pooled.prepare(asas.read_table())
    return output("asas-pooled", {"points": int(time.size), "shift_tests": pooled.subset_tests(pooled.Light(time, flux, pooled.ephemerides(events)))})


@command
def neowise(arguments):
    """The NEOWISE-R L1b W1 and W2 exposures of the star: the visit that holds the eclipse that only 68.95 d predicts (the time of its faintest W2 exposure against the prediction, the window exposures against the others) and the shifts of the comparison stars of the same frames; write out/neowise.json."""
    from hd170312 import neowise as module

    rows = module.read_exposures()
    shifts = module.comparison_shifts(rows)
    return output("neowise", {**module.visit_numbers(rows), "comparison_stars": len(shifts), "comparison_largest_shift": max(shifts.values()), "comparison_shifts": shifts, "visit_dips": module.dip_null(rows)})


ORDER = ["fetch", "eclipses", "radius", "flagged", "qlp", "search", "alias", "forecast", "asas-preregistered", "asas-pooled", "asassn", "neowise", "figure"]


@command
def all(arguments):
    """Run every command that reproduces a claim, in order, then verify the outputs against claims.toml."""
    for name in ORDER:
        print(f"== {name}", file=sys.stderr)
        arguments.record = False
        if COMMANDS[name](arguments):
            return 1
    return verify(arguments)


@command
def verify(arguments):
    """Compare the outputs in out/ with the claims of claims.toml; exit status 1 when a claim fails or its output is missing."""
    from hd170312 import claims

    loaded = claims.load()
    results = claims.verify(loaded) + [("data", *claims.data_check())]
    for line in claims.report(loaded, results[:-1]):
        print(line)
    state, text = results[-1][1:]
    print(f"\nThe outputs were computed from the data files that data/checksums.json records: {state} ({text})")
    failed = [name for name, status, _ in results if status in ("FAIL", "MISSING")]
    skipped = sum(status == "SKIPPED" for _, status, _ in results)
    held = len(results) - len(failed) - skipped
    print(f"{held} of {len(results)} claims hold" + (f", {skipped} skipped" if skipped else "") + (f", {len(failed)} failed" if failed else ""))
    return 1 if failed else 0


def output(name, result):
    """Write a result as out/NAME.json and print it."""
    import json
    from pathlib import Path

    from hd170312 import data

    path = Path(__file__).resolve().parents[2] / "out" / f"{name}.json"
    path.parent.mkdir(exist_ok=True)
    text = json.dumps({**result, "data_sha256": data.data_hashes()} if isinstance(result, dict) else result, indent=1, default=float)
    path.write_text(text + "\n")
    print(text)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="hd170312", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, function in COMMANDS.items():
        p = sub.add_parser(name, help=function.__doc__.split(". ")[0], description=function.__doc__)
        if name == "fetch":
            p.add_argument("--record", action="store_true", help="record the checksum of a file that has none")
    arguments = parser.parse_args(sys.argv[1:] if argv is None else argv)
    return COMMANDS[arguments.command](arguments)
