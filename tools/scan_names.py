"""Scan the tracked files, the staged files and every commit message for internal names (project, experiment, task and role names). The history is public, so a name that enters it stays.

usage: python tools/scan_names.py            every tracked file and every commit message
       python tools/scan_names.py --staged   the staged files (the pre-commit hook)
       python tools/scan_names.py --message FILE   a commit message (the commit-msg hook)"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# Patterns are written in pieces so that this file does not contain the names it forbids.
PATTERNS = {
    "a project name": "star" + r"lord",
    "an experiment number": r"\b[Ee]xperiments?[ /]\d{2,4}\b|\bexperiments/",
    "a task number": r"\bqueue item\b|\bW\d{3}\b",
    "a role name": r"\b(?:agents?|subagents?|manager|referee)\b",
    "a decision attributed to a person in a role": r"\bthe owner\b",
    "the multiplication sign": "×",
}
SKIP_SUFFIXES = {".png", ".jpg", ".jpeg", ".pdf", ".npz", ".fits", ".gz", ".lock"}
SKIP_FILES = {"tools/scan_names.py", "tests/test_names.py"}  # they hold the patterns


def found_in(text):
    return [(what, match.group(0)) for what, pattern in PATTERNS.items() for match in re.finditer(pattern, text)]


def git(*arguments):
    return subprocess.run(["git", *arguments], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def has_commits():
    return subprocess.run(["git", "rev-parse", "--verify", "-q", "HEAD"], cwd=ROOT, capture_output=True, check=False).returncode == 0


def tracked(staged=False):
    names = git("diff", "--cached", "--name-only", "--diff-filter=ACMR") if staged else git("ls-files")
    return [name for name in names.splitlines() if Path(name).suffix not in SKIP_SUFFIXES and name not in SKIP_FILES]


def file_text(name, staged=False):
    if staged:
        return subprocess.run(["git", "show", f":{name}"], cwd=ROOT, capture_output=True, text=True, check=False).stdout
    path = ROOT / name
    return path.read_text(errors="replace") if path.exists() else ""


def problems(staged=False, messages=True):
    found = [(name, what, text) for name in tracked(staged) for what, text in found_in(file_text(name, staged))]
    if messages and not staged and has_commits():
        for commit in git("log", "--format=%H").split():
            found += [(f"commit {commit[:10]}", what, text) for what, text in found_in(git("log", "-1", "--format=%B", commit))]
    return found


def main(argv):
    if argv[:1] == ["--message"]:
        found = found_in(Path(argv[1]).read_text())
        where = "the commit message"
        found = [(where, what, text) for what, text in found]
    else:
        found = problems(staged=argv[:1] == ["--staged"])
    for where, what, text in found:
        print(f"{where}: {what}: {text!r}")
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
