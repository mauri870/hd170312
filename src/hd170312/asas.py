"""The ASAS-3 light curve of 183419-6849.9 (Pojmanski 1997, Acta Astronomica 47, 467; https://www.astrouw.edu.pl/asas/): the archive file and its nightly means.

The file has one block for each field dataset, with columns HJD-2450000, MAG_3, MAG_0, MAG_1, MAG_2, MAG_4, MER_*, GRADE, FRAME; the order of the apertures differs between datasets, so columns are read by name. A magnitude
of 29.999 is not measured. A nightly mean is the mean time and magnitude of the points of one integer HJD day; its error is the robust scatter of the dataset's nightly means over the root of the points of the night.
`read_rows` takes one aperture and the grades asked for (MAG_3, grade A for the pre-registered test); `pooled.prepare` takes grades A and B in the aperture of lowest median frame error."""

import numpy as np

from hd170312.data import ASAS_DESIGNATION, DATA, verified

BTJD_OFFSET = 7000.0  # HJD - 2450000 to HJD - 2457000, the offset of BTJD
NOT_MEASURED = 29.9


def read_table(path=None):
    """The rows of the star: [{hjd (minus 2450000), mag {aperture: value}, err {aperture: value}, grade, dataset}], the columns taken from the header line of each dataset."""
    path = verified(path) if path else verified(DATA / "asas" / f"{ASAS_DESIGNATION}.txt")
    rows, designation, dataset, names = [], None, None, None
    for line in path.read_text().splitlines():
        if line.startswith("#dataset="):
            dataset = line.split("=")[1].strip()
        elif line.startswith("#desig="):
            designation = line.split("=")[1].strip()
        elif line.startswith("#") and "HJD" in line and "MAG_" in line:
            names = line.lstrip("#").split()
        elif line and not line.startswith("#") and designation == ASAS_DESIGNATION:
            fields = line.split()
            mags = {n: float(fields[i]) for i, n in enumerate(names) if n.startswith("MAG_")}
            errs = {n.replace("MER", "MAG"): float(fields[i]) for i, n in enumerate(names) if n.startswith("MER_")}
            rows.append({"hjd": float(fields[0]), "mag": mags, "err": errs, "grade": fields[names.index("GRADE")], "dataset": dataset})
    return rows


def read_rows(path=None, grades=("A",), aperture="MAG_3"):
    """{dataset: array of (HJD - 2457000, magnitude)} for the blocks of the star: the rows of the given grades whose magnitude in the aperture (read by its name) is measured; blocks of one dataset are merged."""
    rows = {}
    for row in read_table(path):
        if row["grade"] in grades and row["mag"][aperture] < NOT_MEASURED:
            rows.setdefault(int(row["dataset"].split()[0]), []).append((row["hjd"] - BTJD_OFFSET, row["mag"][aperture]))
    return {dataset: np.array(values) for dataset, values in rows.items()}


def nightly_means(rows):
    """{dataset: {t, m, c, sc, e}}: the nightly means (time, magnitude, points), the robust scatter of the means and the error of each night."""
    nights = {}
    for dataset, values in rows.items():
        night = np.floor(values[:, 0] + BTJD_OFFSET).astype(int)
        t, m, c = [], [], []
        for k in np.unique(night):
            selected = night == k
            t.append(values[selected, 0].mean())
            m.append(values[selected, 1].mean())
            c.append(selected.sum())
        t, m, c = map(np.array, (t, m, c))
        scatter = 1.4826 * np.median(abs(m - np.median(m)))
        nights[dataset] = {"t": t, "m": m, "c": c, "sc": scatter, "e": scatter / np.sqrt(c)}
    return nights


def stacked(nights):
    """The nights of all datasets as arrays (time, magnitude, error, dataset, points), in the order of the datasets."""
    datasets = sorted(nights)
    return tuple(np.concatenate([nights[d][key] for d in datasets]) for key in ("t", "m", "e")) + (np.concatenate([np.full(len(nights[d]["t"]), d) for d in datasets]), np.concatenate([nights[d]["c"] for d in datasets]))
