from pathlib import Path

import numpy as np
import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def windows():
    """A loader of the fixture windows of the real light curves: sector -> {time, pdcsap, sap, background}."""
    arrays = np.load(FIXTURES / "windows.npz")
    return lambda sector: {name: arrays[f"{name}_{sector}"] for name in ("time", "pdcsap", "sap", "background")}
