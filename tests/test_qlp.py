"""The QLP check on synthetic light curves: the file path, the quality cut, the comparison with SPOC fits and the baseline ramp."""
import numpy as np
from astropy.io import fits

from hd170312 import claims, data, qlp


def write_qlp(directory, sector, time, flux, quality, column="SAP_FLUX"):
    path = directory / "qlp" / data.qlp_path(sector).rsplit("/", 1)[1]
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [fits.Column(name="TIME", format="D", array=time), fits.Column(name=column, format="E", array=flux), fits.Column(name="QUALITY", format="J", array=quality)]
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns)]).writeto(path)


def test_the_path_of_a_qlp_file_is_laid_out_by_sector_and_tic():
    assert data.qlp_path(102) == "s0102/0000/0002/9229/2315/hlsp_qlp_tess_ffi_s0102-0000000292292315_tess_v01_llc.fits"


def test_the_qlp_loader_keeps_quality_zero_and_normalises_by_the_median(tmp_path, monkeypatch):
    monkeypatch.setattr(data, "read_checksums", dict)
    time = np.arange(10.0)
    write_qlp(tmp_path, 102, time, np.array([2.0] * 8 + [50.0, np.nan]), np.array([0] * 8 + [0, 0]))
    write_qlp(tmp_path / "other", 93, time, np.full(10, 3.0), np.array([0, 4096] + [0] * 8))
    columns = data.load_qlp(102, "SAP_FLUX", tmp_path)
    assert columns["time"].size == 9 and abs(np.median(columns["flux"]) - 1.0) < 1e-12
    assert data.load_qlp(93, "SAP_FLUX", tmp_path / "other")["time"].size == 9


def test_the_baseline_ramp_is_the_rise_of_the_medians_across_the_eclipse_and_zero_for_a_flat_flux():
    time = np.arange(0.0, 6.0, 1 / 72)
    flat = qlp.baseline(time, np.ones_like(time), 3.0)
    ramp = qlp.baseline(time, 1.0 + 0.002 * (time - 3.0), 3.0)
    assert abs(flat["ramp"]) < 1e-12 and flat["low"] == flat["high"] == 1.0
    assert 0.002 < ramp["ramp"] < 0.004 and ramp["low"] < 1.0 < ramp["high"]


def test_an_informational_claim_shows_its_value_and_never_fails(tmp_path):
    import json

    (tmp_path / "x.json").write_text(json.dumps({"s": "text"}))
    assert claims.verify([{"id": "i", "command": "x", "path": ["s"], "informational": True}], tmp_path) == [("i", "INFO", "text")]


def test_the_qlp_time_is_compared_in_sector_39_alone_because_the_times_of_sectors_99_to_104_carry_a_barycentric_error(monkeypatch):
    spoc = {"S39": {"time": 2377.1543, "error": 0.00015, "depth": 0.085}, "S102": {"time": 4131.0663, "error": 0.0001, "depth": 0.154}, "S103": {"time": 4169.8482, "error": 0.0005, "depth": 0.085}}
    qlp_times = {"S39": 2377.1546, "S102": 4131.0700, "S103": 4169.8511}  # the last two are shifted by minutes, as the affected sectors are
    monkeypatch.setattr(qlp, "fit_events", lambda directory: spoc)
    monkeypatch.setattr(data, "load_qlp", lambda sector, flux="SAP_FLUX", directory=None: {"time": np.arange(0.0, 10.0, 0.1), "flux": np.ones(100)})
    monkeypatch.setattr(data, "load_sector", lambda sector, directory=None: {"time": np.arange(0.0, 10.0, 0.1), "pdcsap": np.ones(100), "sap": np.ones(100)})
    sectors = iter(["S39", "S102", "S103"])

    def fit(time, flux, guess):
        name = next(sectors)
        return {"time": qlp_times[name], "error": 0.0002, "depth": spoc[name]["depth"]}

    monkeypatch.setattr(qlp, "fit_event_time", fit)
    result = qlp.run()
    assert result["sector_39_time_difference_over_error"] == np.float64(abs(2377.1546 - 2377.1543) / np.hypot(0.0002, 0.00015))
    assert "largest_time_difference_minutes" not in result and "barycentric" in result["summary"], "the other sectors' times are not compared, and the summary says why"
    assert result["comparison"]["S103"]["time_difference_over_error"] > 3, "they are still in the comparison table"
