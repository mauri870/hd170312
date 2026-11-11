"""The public data: the SPOC 2-minute and QLP light curves of TIC 292292315 and the ASAS-3 light curve, fetched, checksummed and read.

SPOC files are in the public bucket `stpubdata` (https://registry.opendata.aws/tess/), listed with S3 ListObjectsV2; QLP files come from the MAST download service (https://archive.stsci.edu/hlsp/qlp). A file is never
replaced; data/checksums.json holds the sha256 of each and a file that does not match is refused."""
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from xml.etree import ElementTree

import numpy as np

TIC = 292292315
SECTORS = (13, 39, 66, 93, 94, 100, 101, 102, 103)
QLP_SECTORS = (39, 93, 102, 103)
ASAS_DESIGNATION = "183419-6849.9"
BUCKET = "https://stpubdata.s3.amazonaws.com"
MAST_DOWNLOAD = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/qlp/"
ASAS_URL = "https://www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data?{designation},asas3"
ASASSN_FILE = "asassn/asassn_v1_saturated.csv"  # the Sky Patrol v1 photometry of saturated stars, downloaded by hand (see `ensure_manual`)
IRSA_TAP = "https://irsa.ipac.caltech.edu/TAP/sync"  # https://irsa.ipac.caltech.edu/docs/program_interface/TAP.html
NEOWISE_TABLE = "neowiser_p1bs_psd"  # the NEOWISE-R Single Exposure (L1b) Source Table
NEOWISE_COLUMNS = "ra,dec,mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro,cc_flags,ph_qual,qual_frame"
NEOWISE_FIELD_COLUMNS = "ra,dec,mjd,w1mpro,w2mpro,cc_flags,ph_qual,qual_frame,allwise_cntr,w2mpro_allwise"
STAR = (278.57828, -68.83253)  # RA, Dec (deg) of HD 170312
FIELD_MJD = (56922, 56927)  # the days of the NEOWISE visit with the eclipse that only 68.95 d predicts
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
CHECKSUMS = DATA / "checksums.json"
_NAMESPACE = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}


class ChecksumMismatch(RuntimeError):
    """A file does not match the sha256 recorded for it."""


def prefix(tic, sector):
    """The bucket prefix of the products of a star in a sector."""
    digits = f"{int(tic):016d}"
    return f"tess/public/tid/s{sector:04d}/{digits[0:4]}/{digits[4:8]}/{digits[8:12]}/{digits[12:16]}/"


def parse_listing(content):
    """(keys, continuation token or None) of one page of a ListObjectsV2 response."""
    root = ElementTree.fromstring(content)
    keys = [element.text for element in root.findall("s3:Contents/s3:Key", _NAMESPACE) if element.text]
    token = root.find("s3:NextContinuationToken", _NAMESPACE)
    return keys, token.text if token is not None and token.text else None


def http_get(url, attempts=4, timeout=300):
    """The bytes of a URL, retried with a growing pause."""
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read()
        except OSError:
            if attempt == attempts - 1:
                raise
            time.sleep(2 ** attempt)


def light_curve_key(sector, get=http_get):
    """The bucket key of the SPOC 2-minute light curve of the star in a sector."""
    keys, token = [], None
    while True:
        query = {"list-type": "2", "prefix": prefix(TIC, sector)}
        if token:
            query["continuation-token"] = token
        page, token = parse_listing(get(f"{BUCKET}/?{urllib.parse.urlencode(query)}"))
        keys += page
        if token is None:
            break
    found = [key for key in keys if key.endswith("-s_lc.fits")]
    if len(found) != 1:
        raise LookupError(f"sector {sector}: {len(found)} light curve files in {prefix(TIC, sector)}")
    return found[0]


def qlp_path(sector):
    """The MAST path of the QLP light curve of the star in a sector, laid out by sector and TIC identifier (https://archive.stsci.edu/hlsp/qlp)."""
    digits = f"{TIC:016d}"
    return f"s{sector:04d}/{digits[0:4]}/{digits[4:8]}/{digits[8:12]}/{digits[12:16]}/hlsp_qlp_tess_ffi_s{sector:04d}-{digits}_tess_v01_llc.fits"


def neowise_url(columns, radius_arcsec, where=""):
    """The IRSA TAP query of the rows of the L1b table within `radius_arcsec` of the star, ordered by time and by position, as CSV."""
    cone = f"contains(point('ICRS',ra,dec),circle('ICRS',{STAR[0]},{STAR[1]},{radius_arcsec / 3600:.7f}))=1"
    query = f"select {columns} from {NEOWISE_TABLE} where {cone}{where} order by mjd, ra"
    return f"{IRSA_TAP}?{urllib.parse.urlencode({'QUERY': query, 'FORMAT': 'csv'})}"


def neowise_target_url():
    return neowise_url(NEOWISE_COLUMNS, 10)


def neowise_field_url():
    return neowise_url(NEOWISE_FIELD_COLUMNS, 600, f" and mjd between {FIELD_MJD[0]} and {FIELD_MJD[1]}")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_checksums(path=CHECKSUMS):
    return json.loads(Path(path).read_text()) if Path(path).exists() else {}


def check(path, name, checksums):
    """Raise ChecksumMismatch when the file differs from its recorded sha256; return whether a checksum was recorded for it."""
    recorded = checksums.get(name)
    if recorded is None:
        return False
    found = sha256(path)
    if found != recorded["sha256"]:
        raise ChecksumMismatch(f"{name}: sha256 {found}, recorded {recorded['sha256']}")
    return True


def ensure(name, url, directory=DATA, checksums=None, record=False, get=http_get):
    """The path of a data file, downloaded when absent and checked against data/checksums.json; with `record` a file with no recorded checksum is added to it."""
    checksums = read_checksums() if checksums is None else checksums
    path = Path(directory) / name
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        content = get(url)
        path.write_bytes(content)
    known = check(path, name, checksums)
    if not known and record:
        checksums[name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    return path


def ensure_manual(name, directory=DATA, checksums=None):
    """The path of a data file that cannot be fetched by a script, checked against data/checksums.json; a file that is not there is an error that says how to get it.

    The ASAS-SN photometry of saturated stars (Winecki & Kochanek 2024; Shappee et al. 2014, ApJ 788, 48; Kochanek et al. 2017, PASP 129, 104502) comes from the form "Compute Coordinate Lightcurve" of
    https://asas-sn.osu.edu/ (Sky Patrol v1) with the method "Saturated Stars (Machine Learning)", RA 278.57828, Dec -68.83253, which asks for a reCAPTCHA. The file in the repository was downloaded on 2026 November 11."""
    checksums = read_checksums() if checksums is None else checksums
    path = Path(directory) / name
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing: it is tracked in the repository; to make it again, submit the form 'Compute Coordinate Lightcurve' of https://asas-sn.osu.edu/ (Saturated Stars (Machine Learning)) at RA {STAR[0]}, Dec {STAR[1]}")
    check(path, name, checksums)
    return path


def fetch(record=False, directory=DATA, get=http_get):
    """Download the nine SPOC light curves, the four QLP light curves, the ASAS-3 file and the two NEOWISE queries, and check the ASAS-SN file that the repository holds; return {name: path}."""
    checksums = read_checksums()
    paths = {}
    for sector in SECTORS:
        key = light_curve_key(sector, get)
        name = f"tess/{Path(key).name}"
        paths[name] = ensure(name, f"{BUCKET}/{key}", directory, checksums, record, get)
    for sector in QLP_SECTORS:
        relative = qlp_path(sector)
        name = f"qlp/{Path(relative).name}"
        paths[name] = ensure(name, MAST_DOWNLOAD + relative, directory, checksums, record, get)
    name = f"asas/{ASAS_DESIGNATION}.txt"
    paths[name] = ensure(name, ASAS_URL.format(designation=ASAS_DESIGNATION), directory, checksums, record, get)
    paths[ASASSN_FILE] = ensure_manual(ASASSN_FILE, directory, checksums)
    for name, url in (("neowise/l1b_10arcsec.csv", neowise_target_url()), ("neowise/l1b_visit_field_600arcsec.csv", neowise_field_url())):
        paths[name] = ensure(name, url, directory, checksums, record, get)
    if record:
        CHECKSUMS.write_text(json.dumps(dict(sorted(checksums.items())), indent=1) + "\n")
    return paths


def verified(path, directory=DATA, checksums=None):
    """The path of a data file after its sha256 is checked against the recorded one (a file under `directory` with no recorded checksum is returned unchecked, as in the tests, whose files are synthetic)."""
    checksums = read_checksums() if checksums is None else checksums
    path = Path(path)
    try:
        name = path.resolve().relative_to(Path(directory).resolve()).as_posix()
    except ValueError:
        return path
    check(path, name, checksums)
    return path


def data_hashes(directory=DATA):
    """{name: sha256} of the data files that are recorded in data/checksums.json and present, each checked against the record; the data an output came from."""
    checksums = read_checksums()
    hashes = {}
    for name in sorted(checksums):
        path = Path(directory) / name
        if path.exists():
            check(path, name, checksums)
            hashes[name] = checksums[name]["sha256"]
    return hashes


def tess_file(sector, directory=DATA):
    """The path of the light curve file of a sector that `fetch` downloaded, checked against its recorded sha256."""
    found = sorted((Path(directory) / "tess").glob(f"tess*-s{sector:04d}-{TIC:016d}-*-s_lc.fits"))
    if len(found) != 1:
        raise FileNotFoundError(f"no light curve of sector {sector} in {Path(directory) / 'tess'}: run `hd170312 fetch`")
    return verified(found[0], directory)


def load_sector(sector, directory=DATA):
    """The columns of a sector's light curve, cadences of QUALITY 0 with a finite time and PDCSAP flux: time (BTJD), pdcsap (relative to the sector median), sap, background (e/s) and the CROWDSAP keyword."""
    from astropy.io import fits

    with fits.open(tess_file(sector, directory)) as hdus:
        table = hdus[1].data
        good = (table["QUALITY"] == 0) & np.isfinite(table["TIME"]) & np.isfinite(table["PDCSAP_FLUX"])
        pdcsap = np.asarray(table["PDCSAP_FLUX"][good], dtype=float)
        return {
            "time": np.asarray(table["TIME"][good], dtype=float),
            "pdcsap": pdcsap / np.median(pdcsap),
            "sap": np.asarray(table["SAP_FLUX"][good], dtype=float),
            "background": np.asarray(table["SAP_BKG"][good], dtype=float),
            "crowdsap": float(hdus[1].header.get("CROWDSAP", np.nan)),
        }


def load_sector_sap(sector, directory=DATA):
    """The SAP flux of a sector with the cadences flagged for scattered light (QUALITY bit 4096) kept where the flux is finite: PDCSAP is empty there, and the eclipse of sector 93 and the egress of sector 101 lie in those cadences; time (BTJD), sap (e/s), background (e/s) and quality."""
    from astropy.io import fits

    with fits.open(tess_file(sector, directory)) as hdus:
        table = hdus[1].data
        good = np.isin(table["QUALITY"], (0, 4096)) & np.isfinite(table["TIME"]) & np.isfinite(table["SAP_FLUX"])
        return {name: np.asarray(table[column][good], dtype=float if name != "quality" else int) for name, column in (("time", "TIME"), ("sap", "SAP_FLUX"), ("background", "SAP_BKG"), ("quality", "QUALITY"))}


def load_qlp(sector, flux="SAP_FLUX", directory=DATA):
    """The QLP light curve of a sector: time (BTJD) and the flux column divided by its median, for the cadences of QUALITY 0 with a finite time and flux."""
    from astropy.io import fits

    path = verified(Path(directory) / "qlp" / Path(qlp_path(sector)).name, directory)
    with fits.open(path) as hdus:
        table = hdus[1].data
        good = (table["QUALITY"] == 0) & np.isfinite(table["TIME"]) & np.isfinite(table[flux])
        values = np.asarray(table[flux][good], dtype=float)
        return {"time": np.asarray(table["TIME"][good], dtype=float), "flux": values / np.median(values)}
