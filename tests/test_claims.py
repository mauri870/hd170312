"""claims.toml names commands and values that exist, and `verify` compares outputs with it."""
import json

from hd170312 import claims
from hd170312.cli import COMMANDS, ORDER

CLAIMS = claims.load()


def test_every_claim_names_a_command_of_the_run_and_has_a_value_to_compare():
    for claim in CLAIMS:
        assert claim["command"] in COMMANDS and claim["command"] in ORDER and claim["path"], claim["id"]
        assert claim.get("informational") or "expected" in claim and (claim.get("kind") in ("equals", "at_most") or "tolerance" in claim), claim["id"]


def test_the_ids_of_the_claims_are_unique():
    ids = [claim["id"] for claim in CLAIMS]
    assert len(ids) == len(set(ids))


def test_a_value_is_within_its_tolerance_after_the_scale_and_a_list_or_a_string_is_equal():
    value = {"expected": 1.0, "tolerance": 0.1}
    assert claims.check(value, 1.05)[0] and not claims.check(value, 1.2)[0]
    assert claims.check({"expected": 8.5, "tolerance": 0.05, "scale": 100.0}, 0.0852)[0]
    assert claims.check({"kind": "equals", "expected": [1, 2]}, [1, 2])[0] and not claims.check({"kind": "equals", "expected": [1, 2]}, [1])[0]


def test_verify_reports_pass_fail_missing_and_info(tmp_path):
    (tmp_path / "x.json").write_text(json.dumps({"a": {"b": [3.0, 4.0]}, "s": "text"}))
    rows = [{"id": "ok", "command": "x", "path": ["a", "b", "1"], "expected": 4.0, "tolerance": 0.1}, {"id": "bad", "command": "x", "path": ["a", "b", "0"], "expected": 4.0, "tolerance": 0.1},
            {"id": "gone", "command": "y", "path": ["a"], "expected": 1.0, "tolerance": 0.1}, {"id": "bad-path", "command": "x", "path": ["zzz"], "expected": 1.0, "tolerance": 0.1},
            {"id": "info", "command": "x", "path": ["s"], "informational": True}]
    status = {name: s for name, s, _ in claims.verify(rows, tmp_path)}
    assert status == {"ok": "PASS", "bad": "FAIL", "gone": "MISSING", "bad-path": "FAIL", "info": "INFO"}


def test_the_data_check_passes_matching_hashes_and_fails_other_data_and_a_missing_record(tmp_path, monkeypatch):
    from hd170312 import data

    monkeypatch.setattr(data, "read_checksums", lambda: {"a.txt": {"sha256": "aa", "bytes": 1}})
    (tmp_path / "ok.json").write_text(json.dumps({"data_sha256": {"a.txt": "aa"}}))
    assert claims.data_check(tmp_path)[0] == "PASS"
    (tmp_path / "other.json").write_text(json.dumps({"data_sha256": {"a.txt": "bb"}}))
    (tmp_path / "none.json").write_text(json.dumps({"x": 1}))
    status, text = claims.data_check(tmp_path)
    assert status == "FAIL" and "other.json" in text and "none.json" in text and "ok.json" not in text


def test_without_the_group_the_search_is_skipped_and_its_claims_are_neither_passed_nor_failed(monkeypatch, tmp_path, capsys):
    from hd170312 import cli

    monkeypatch.setattr(claims, "search_installed", lambda: False)
    assert cli.search(None) == 0 and "--group search" in capsys.readouterr().err
    rows = [{"id": "p", "command": "search", "requires": "search", "path": ["sde"], "expected": 10.0, "tolerance": 1.0}]
    assert claims.verify(rows, tmp_path) == [("p", "SKIPPED", claims.SEARCH_HINT)]


def test_every_claim_has_a_statement_and_one_of_the_questions():
    assert len(claims.QUESTIONS) == 11 and len(set(claims.QUESTIONS)) == 11
    for claim in CLAIMS:
        assert claim.get("statement", "").strip() and not claim["statement"].endswith("."), claim["id"]
        assert claim.get("question") in claims.QUESTIONS, claim["id"]
    assert {claim["question"] for claim in CLAIMS} == set(claims.QUESTIONS), "every question has a claim under it"


def test_the_report_groups_the_claims_under_their_questions_in_the_order_of_the_questions():
    rows = [{"id": "b", "question": claims.QUESTIONS[2], "statement": "Second group"}, {"id": "a", "question": claims.QUESTIONS[0], "statement": "First group"},
            {"id": "c", "question": claims.QUESTIONS[0], "statement": "First group again"}]
    results = [("b", "PASS", "2"), ("a", "PASS", "1"), ("c", "PASS", "3")]
    lines = claims.report(rows, results)
    assert [line for line in lines if line.endswith("?")] == [claims.QUESTIONS[0], claims.QUESTIONS[2]]
    assert [line.split()[1] for line in lines if line.startswith("  ")] == ["First", "First", "Second"], "the order of the file inside a group, the groups in the order of the questions"


def test_a_line_gives_the_status_the_statement_and_the_computed_value_and_the_tolerance_only_on_failure():
    rows = [{"id": "ok", "question": claims.QUESTIONS[0], "statement": "A is 4"}, {"id": "bad", "question": claims.QUESTIONS[0], "statement": "B is 4"},
            {"id": "gone", "question": claims.QUESTIONS[1], "statement": "C is 1"}, {"id": "skip", "question": claims.QUESTIONS[5], "statement": "D is 2"}]
    results = [("ok", "PASS", "4"), ("bad", "FAIL", "3 against 4 +- 0.1"), ("gone", "MISSING", "run `hd170312 y`"), ("skip", "SKIPPED", "uv sync --locked --group search")]
    lines = [line.strip() for line in claims.report(rows, results) if line.startswith("  ")]
    assert lines == ["PASS     A is 4 (computed 4)", "FAIL     B is 4 (computed 3 against 4 +- 0.1)", "MISSING  C is 1 (run `hd170312 y`)", "SKIPPED  D is 2 (uv sync --locked --group search)"]


def test_check_gives_the_computed_value_apart_from_the_comparison():
    assert claims.check({"expected": 8.5, "tolerance": 0.05, "scale": 100.0}, 0.0852) == (True, "8.52", "8.52 against 8.5 +- 0.05")
    assert claims.check({"kind": "at_most", "expected": 1.6}, 2.0)[:2] == (False, "2.00")
    assert claims.check({"kind": "equals", "expected": "2027-03-11"}, "2027-03-12")[1] == "'2027-03-12'"


def test_the_computed_value_is_rounded_to_the_place_of_the_tolerance_and_integers_and_strings_stay_exact():
    computed = lambda claim, value: claims.check(claim, value)[1]
    assert computed({"expected": 2377.1543, "tolerance": 5e-05}, 2377.154261) == "2377.15426"
    assert computed({"expected": 8.5, "tolerance": 0.05, "scale": 100.0}, 0.08546003277) == "8.55"
    assert computed({"expected": 1.03, "tolerance": 0.005}, 1.025904816) == "1.026"
    assert computed({"expected": 0.011, "tolerance": 0.0005}, 0.011099778) == "0.0111"
    assert computed({"expected": 0.0049, "tolerance": 0.00005}, 0.004901670202) == "0.00490"
    assert computed({"expected": 3193, "tolerance": 0}, 3193) == "3193", "an exact integer stays"
    assert computed({"kind": "equals", "expected": "2027-03-11"}, "2027-03-11") == "'2027-03-11'"
    assert computed({"kind": "at_most", "expected": 1.6}, 1.577296409) == "1.58", "an at_most claim shows one more digit than its limit"
    assert computed({"kind": "at_most", "expected": 0.001}, 0.0002281466495) == "0.0002"
    assert computed({"kind": "at_most", "expected": 4.1}, 4.053005795) == "4.05"
    assert claims.check({"expected": 8.5, "tolerance": 0.05, "scale": 100.0}, 0.0852)[2] == "8.52 against 8.5 +- 0.05", "the comparison of a failing claim keeps its value"


def test_the_neowise_comparison_claim_limits_the_largest_fading_to_the_0_02_mag_that_the_note_states():
    claim = next(c for c in CLAIMS if c["id"] == "neowise-comparison-largest-shift")
    assert claim["kind"] == "at_most" and claim["expected"] == 0.02 and "0.02" in claim["statement"]
