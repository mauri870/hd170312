"""Write tests/fixtures/windows.npz from the fetched light curves: the PDCSAP flux, the SAP flux and the background in 6 d around each of the three complete eclipses (sectors 39, 102, 103) and in the last 3 d of sector 93, which the offline tests use in place of the full files.

usage: uv run python tools/make_fixtures.py   (after `uv run hd170312 fetch`)"""
from pathlib import Path

import numpy as np

from hd170312.data import load_sector

OUT = Path(__file__).resolve().parent.parent / "tests" / "fixtures"
WINDOWS = {39: (2374.16, 2380.16), 102: (4128.07, 4134.07), 103: (4166.85, 4172.85), 93: (3849.5, 3855.1)}  # BTJD


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    arrays = {}
    for sector, (low, high) in WINDOWS.items():
        columns = load_sector(sector)
        inside = (columns["time"] >= low) & (columns["time"] <= high)
        for name in ("time", "pdcsap", "sap", "background"):
            arrays[f"{name}_{sector}"] = columns[name][inside].astype(np.float64)
        arrays[f"crowdsap_{sector}"] = np.array(columns["crowdsap"])
    np.savez_compressed(OUT / "windows.npz", **arrays)
    print(f"wrote {OUT / 'windows.npz'}: {(OUT / 'windows.npz').stat().st_size} bytes")


if __name__ == "__main__":
    main()
