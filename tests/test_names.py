"""No tracked file and no commit message names a project, an experiment, a task or a role (tools/scan_names.py); the check runs before each commit (tools/install_hooks.sh) and here."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import scan_names


@pytest.mark.parametrize("text, what", [
    ("ported from " + "star" + "lord", "a project name"),
    ("see experiment 123 and experiments/456", "an experiment number"),
    ("queue item W356", "a task number"),
    ("the manager and a fresh agent read it", "a role name"),
    ("as the owner decided", "a decision attributed to a person in a role"),
    ("3 × 4", "the multiplication sign"),
])
def test_every_kind_of_name_is_found(text, what):
    assert what in [kind for kind, _ in scan_names.found_in(text)]


def test_numbers_and_public_names_are_not_found():
    assert scan_names.found_in("A = 1.035 +- 0.048 from TIC 292292315 with MAST, Gaia DR3 and ASAS-3; 2 x 10^-5; the experiment ran once") == []


def test_no_tracked_file_and_no_commit_message_names_it():
    assert scan_names.problems() == []
