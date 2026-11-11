"""The data layer: bucket layout, listing, checksums and the reading of a light curve file, offline."""
import json

import numpy as np
import pytest
from astropy.io import fits

from hd170312 import data

LISTING = """<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/"><Name>stpubdata</Name><IsTruncated>false</IsTruncated>
<Contents><Key>tess/public/tid/s0039/0000/0002/9229/2315/tess2021146024351-s0039-0000000292292315-0210-s_lc.fits</Key></Contents>
<Contents><Key>tess/public/tid/s0039/0000/0002/9229/2315/tess2021146024351-s0039-0000000292292315-0210-s_tp.fits</Key></Contents></ListBucketResult>"""


def test_the_prefix_groups_the_sixteen_digits_of_the_tic_in_four():
    assert data.prefix(292292315, 39) == "tess/public/tid/s0039/0000/0002/9229/2315/"


def test_the_light_curve_key_is_the_one_that_ends_in_s_lc():
    seen = []

    def get(url):
        seen.append(url)
        return LISTING.encode()

    assert data.light_curve_key(39, get).endswith("-s0039-0000000292292315-0210-s_lc.fits")
    assert "prefix=tess%2Fpublic%2Ftid%2Fs0039%2F0000%2F0002%2F9229%2F2315%2F" in seen[0]


def test_a_sector_with_no_light_curve_is_an_error_and_never_a_neighbours_file():
    empty = '<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/"><IsTruncated>false</IsTruncated></ListBucketResult>'
    with pytest.raises(LookupError):
        data.light_curve_key(39, lambda url: empty.encode())


def test_a_file_that_differs_from_its_recorded_checksum_is_refused(tmp_path):
    (tmp_path / "a.txt").write_text("first")
    checksums = {"a.txt": {"sha256": data.sha256(tmp_path / "a.txt"), "bytes": 5}}
    assert data.ensure("a.txt", "http://unused", tmp_path, checksums) == tmp_path / "a.txt"
    (tmp_path / "a.txt").write_text("second")
    with pytest.raises(data.ChecksumMismatch):
        data.ensure("a.txt", "http://unused", tmp_path, checksums)


def test_a_missing_file_is_downloaded_and_a_file_with_no_checksum_is_recorded_on_request(tmp_path):
    checksums = {}
    path = data.ensure("sub/b.txt", "http://unused", tmp_path, checksums, record=True, get=lambda url: b"content")
    assert path.read_bytes() == b"content" and checksums["sub/b.txt"]["bytes"] == 7
    assert data.ensure("sub/b.txt", "http://unused", tmp_path, checksums, get=lambda url: b"other") == path  # an existing file is not fetched again


def test_the_tracked_checksums_name_the_nine_light_curves_and_the_asas_file():
    recorded = json.loads(data.CHECKSUMS.read_text())
    assert sorted(name for name in recorded if name.startswith("tess/")) and len([n for n in recorded if n.startswith("tess/")]) == len(data.SECTORS)
    assert f"asas/{data.ASAS_DESIGNATION}.txt" in recorded and all(len(entry["sha256"]) == 64 for entry in recorded.values())


def test_load_sector_keeps_cadences_of_quality_zero_with_finite_flux_and_normalizes_by_the_median(tmp_path):
    rows = 6
    time = np.array([1.0, 2.0, 3.0, 4.0, np.nan, 6.0])
    pdcsap = np.array([10.0, 20.0, 30.0, np.nan, 50.0, 60.0])
    quality = np.array([0, 0, 0, 0, 0, 8])
    columns = [fits.Column(name=n, format="D", array=a) for n, a in (("TIME", time), ("SAP_FLUX", pdcsap), ("PDCSAP_FLUX", pdcsap), ("SAP_BKG", np.ones(rows)))] + [fits.Column(name="QUALITY", format="J", array=quality)]
    table = fits.BinTableHDU.from_columns(columns)
    table.header["CROWDSAP"] = 0.99
    (tmp_path / "tess").mkdir()
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(tmp_path / "tess" / f"tess2000000000000-s0039-{data.TIC:016d}-0000-s_lc.fits")
    sector = data.load_sector(39, tmp_path)
    assert list(sector["time"]) == [1.0, 2.0, 3.0] and list(sector["pdcsap"]) == [0.5, 1.0, 1.5] and sector["crowdsap"] == 0.99


def test_a_data_file_that_differs_from_its_recorded_checksum_is_refused_by_the_readers(tmp_path):
    (tmp_path / "tess").mkdir()
    name = f"tess/tess2000000000000-s0039-{data.TIC:016d}-0000-s_lc.fits"
    (tmp_path / name).write_bytes(b"not the file")
    checksums = {name: {"sha256": "0" * 64, "bytes": 12}}
    with pytest.raises(data.ChecksumMismatch):
        data.verified(tmp_path / name, tmp_path, checksums)
    assert data.verified(tmp_path / "elsewhere.txt", tmp_path, checksums) == tmp_path / "elsewhere.txt"


def test_the_hashes_of_the_data_are_those_of_the_present_files_and_a_changed_file_is_refused(tmp_path, monkeypatch):
    (tmp_path / "a.txt").write_text("one")
    digest = data.sha256(tmp_path / "a.txt")
    monkeypatch.setattr(data, "read_checksums", lambda path=data.CHECKSUMS: {"a.txt": {"sha256": digest}, "missing.txt": {"sha256": "0" * 64}})
    assert data.data_hashes(tmp_path) == {"a.txt": digest}
    (tmp_path / "a.txt").write_text("two")
    with pytest.raises(data.ChecksumMismatch):
        data.data_hashes(tmp_path)


def test_load_sector_sap_keeps_the_cadences_flagged_for_scattered_light_where_sap_is_finite(tmp_path):
    time = np.array([1.0, 2.0, 3.0, 4.0, 5.0, np.nan, 7.0])
    sap = np.array([10.0, 20.0, 30.0, 40.0, np.nan, 60.0, 70.0])
    pdcsap = np.array([10.0, 20.0, np.nan, np.nan, np.nan, 60.0, 70.0])
    quality = np.array([0, 4096, 4096, 8, 4096, 0, 128])
    columns = [fits.Column(name=n, format="D", array=a) for n, a in (("TIME", time), ("SAP_FLUX", sap), ("PDCSAP_FLUX", pdcsap), ("SAP_BKG", np.arange(7.0)))] + [fits.Column(name="QUALITY", format="J", array=quality)]
    (tmp_path / "tess").mkdir()
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns)]).writeto(tmp_path / "tess" / f"tess2000000000000-s0093-{data.TIC:016d}-0000-s_lc.fits")
    sector = data.load_sector_sap(93, tmp_path)
    assert list(sector["time"]) == [1.0, 2.0, 3.0] and list(sector["quality"]) == [0, 4096, 4096] and list(sector["sap"]) == [10.0, 20.0, 30.0] and list(sector["background"]) == [0.0, 1.0, 2.0]
    assert list(data.load_sector(93, tmp_path)["time"]) == [1.0], "the PDCSAP reader still keeps QUALITY 0 only"


def test_fetch_downloads_the_two_neowise_queries_and_records_their_checksums(tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setattr(data, "CHECKSUMS", tmp_path / "checksums.json")
    monkeypatch.setattr(data, "read_checksums", lambda path=None: {})
    monkeypatch.setattr(data, "light_curve_key", lambda sector, get=None: f"tess/x/tess2000000000000-s{sector:04d}-{data.TIC:016d}-0000-s_lc.fits")

    def get(url):
        seen[url.split("?")[0]] = seen.get(url.split("?")[0], 0) + 1
        return b"synthetic " + url.encode()

    (tmp_path / "asassn").mkdir()
    (tmp_path / data.ASASSN_FILE).write_text("synthetic, as downloaded by hand\n")
    monkeypatch.setattr(data, "read_checksums", lambda path=None: {data.ASASSN_FILE: {"sha256": data.sha256(tmp_path / data.ASASSN_FILE), "bytes": 33}})
    paths = data.fetch(record=True, directory=tmp_path, get=get)
    assert paths[data.ASASSN_FILE] == tmp_path / data.ASASSN_FILE and not any("asassn" in url for url in seen), "the ASAS-SN file is checked and never downloaded"
    assert {"neowise/l1b_10arcsec.csv", "neowise/l1b_visit_field_600arcsec.csv"} <= set(paths)
    recorded = json.loads((tmp_path / "checksums.json").read_text())
    assert recorded["neowise/l1b_10arcsec.csv"]["sha256"] == data.sha256(tmp_path / "neowise/l1b_10arcsec.csv") and "neowise/l1b_visit_field_600arcsec.csv" in recorded
    text = (tmp_path / "neowise/l1b_visit_field_600arcsec.csv").read_text()
    assert "56922" in text and "56927" in text and "0.1666667" in text, "the field is the 600 arcsec cone on the days of the visit"
