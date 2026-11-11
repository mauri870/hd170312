"""claims.toml lists the main claims of the Note. `verify` compares each with the output of its command in out/."""
import importlib.util
import json
import math
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLAIMS = ROOT / "claims.toml"
OUT = ROOT / "out"


SEARCH_HINT = "the box search needs transitgpu: uv sync --locked --group search"

# The questions that the claims answer, in the order that `verify` prints them; the `question` of each claim is one of these.
QUESTIONS = (
    "When did the eclipses happen, and how deep were they?",
    "Which orbital periods fit the eclipses?",
    "When is the next deciding eclipse?",
    "Did ASAS-3 see the star fainter at those times?",
    "Did ASAS-SN see the star fainter at those times?",
    "Did NEOWISE see an eclipse that only 68.95 days predicts?",
    "Is the fourth eclipse in the flagged cadences of sector 93?",
    "Does sector 101 show the egress that only 68.95 days predicts?",
    "Does the discovery search flag the star?",
    "Is the companion a star?",
    "Does a second reduction of the data (QLP) agree?",
)


def search_installed():
    """Whether transitgpu, which the box search needs, is installed (the group `search`)."""
    return importlib.util.find_spec("transitgpu") is not None


def load(path=CLAIMS):
    return tomllib.loads(Path(path).read_text())["claim"]


def resolve(output, path):
    """The value at a path of keys (an integer key indexes a list) in a command's output."""
    for key in path:
        output = output[int(key)] if isinstance(output, list) else output[str(key)]
    return output


def rounded(claim, value):
    """The value as `verify` shows it: an integer of a claim with a tolerance of zero as it is, otherwise to the decimal place of the claim's tolerance (the place of the first digit that is not zero), or, for a claim with a limit and no tolerance, to one more place than the limit has."""
    if claim.get("tolerance") == 0 and float(value).is_integer():
        return str(int(value))
    if claim.get("tolerance"):
        places = -math.floor(math.log10(claim["tolerance"]))
    else:
        places = len(f"{claim['expected']:.10g}".partition(".")[2]) + 1
    return f"{value:.{max(places, 0)}f}"


def check(claim, measured):
    """(passed, computed, comparison) for a measured value: equal to the expected, at most the expected (kind at_most), or within the tolerance of it."""
    expected = claim["expected"]
    if claim.get("kind") == "equals":
        return measured == expected, f"{measured!r}", f"{measured!r} against {expected!r}"
    value = measured * claim.get("scale", 1.0)
    full = f"{value:.10g}"
    if claim.get("kind") == "at_most":
        return value <= expected, rounded(claim, value), f"{full} against at most {expected:.10g}"
    return abs(value - expected) <= claim["tolerance"], rounded(claim, value), f"{full} against {expected:.10g} +- {claim['tolerance']:.3g}"


def verify(claims, out=OUT):
    """[(id, status, text)]: PASS (the text is the computed value) or FAIL (the computed value against the expected one and the tolerance), MISSING when the output file is absent, SKIPPED for a claim whose `requires` group is not installed, INFO for a claim with `informational = true`, whose value is shown and not judged."""
    results, cache = [], {}
    for claim in claims:
        name = claim["command"]
        if claim.get("requires") == "search" and not search_installed():
            results.append((claim["id"], "SKIPPED", SEARCH_HINT))
            continue
        if name not in cache:
            path = Path(out) / f"{name}.json"
            cache[name] = json.loads(path.read_text()) if path.exists() else None
        if cache[name] is None:
            results.append((claim["id"], "MISSING", f"run `hd170312 {name}`"))
            continue
        try:
            measured = resolve(cache[name], claim["path"])
            passed, computed, comparison = (True, str(measured), "") if claim.get("informational") else check(claim, measured)
        except (KeyError, IndexError, ValueError, TypeError) as error:
            results.append((claim["id"], "FAIL", f"could not be computed: {type(error).__name__}: {error}"))
            continue
        results.append((claim["id"], "INFO" if claim.get("informational") else "PASS" if passed else "FAIL", computed if passed else comparison))
    return results


def report(claims, results):
    """The lines that `verify` prints: each question of QUESTIONS, in its order, with one line under it for each claim of that question (in the order of the file): the status, the statement and, in parentheses, what was computed
    (with the expected value and the tolerance when the claim fails) or why it was not."""
    status = {name: (state, text) for name, state, text in results}
    width = max(len(state) for state, _ in status.values())
    lines = []
    for question in QUESTIONS:
        group = [claim for claim in claims if claim["question"] == question]
        if group:
            lines.append(question)
        for claim in group:
            state, text = status[claim["id"]]
            shown = f"computed {text}" if state in ("PASS", "FAIL", "INFO") and not text.startswith("could not") else text
            lines.append(f"  {state:{width}s}  {claim['statement']} ({shown})")
    return lines


def data_check(out=OUT):
    """(status, text): the data files each output was computed from are the files of data/checksums.json."""
    from hd170312 import data

    recorded = {name: entry["sha256"] for name, entry in data.read_checksums().items()}
    wrong = []
    for path in sorted(Path(out).glob("*.json")):
        used = json.loads(path.read_text()).get("data_sha256")
        if used is None or any(recorded.get(name) != digest for name, digest in used.items()):
            wrong.append(path.name)
    return ("FAIL", f"outputs not from the recorded data: {wrong}; run their commands again") if wrong else ("PASS", "every output used the files of data/checksums.json")
