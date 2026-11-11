"""The Note: its source holds no comments, its figure exists, and it builds."""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

NOTE = Path(__file__).resolve().parent.parent / "note"


def test_the_source_of_the_note_has_no_comment_and_names_its_figure():
    text = (NOTE / "main.tex").read_text()
    assert not re.search(r"(?m)^\s*%|(?<!\\)%", text), "a comment in note/main.tex"
    figures = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", text)
    assert figures == ["figures/events.pdf"] and all((NOTE / name).exists() for name in figures), "the Note's figure is the PDF, a format the AAS accepts (PNG is not among them)"


@pytest.mark.skipif(shutil.which("tectonic") is None, reason="tectonic is not installed")
def test_the_note_builds_into_a_pdf(tmp_path):
    done = subprocess.run(["tectonic", "--outdir", str(tmp_path), "main.tex"], cwd=NOTE, capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr[-2000:]
    assert (tmp_path / "main.pdf").stat().st_size > 50_000
