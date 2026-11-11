"""A synthetic ASAS-3 archive file with the cadence and the scatter of the real one (one point a night on about every third night of 3,200 d, 0.010 mag), for the offline tests: the real file is fetched, never committed."""
import numpy as np

from hd170312 import preregistered

DESIGNATION = "183419-6849.9"


def synthetic_file(path, amplitude=1.0, seed=7, scatter=0.010, other_block=True):
    """Write the file; the eclipses of the Note are injected at the epochs of 68.95 d with `amplitude` times the TESS depth. Returns the path."""
    generator = np.random.default_rng(seed)
    lines = ["# The All Sky Automated Survey Data (synthetic)", "# HJD MAG_3 MAG_0 MAG_1 MAG_2 MAG_4 MER_3 MER_0 MER_1 MER_2 MER_4 GRADE FRAME"]
    for dataset, offset in ((1, 8.71), (2, 8.72), (3, 8.73)):
        lines += [f"#dataset= {dataset} ; 1 F{dataset}", f"#desig= {DESIGNATION}"]
        nights = np.sort(generator.choice(np.arange(-5130, -1895), size=380, replace=False)).astype(float) + generator.uniform(0.1, 0.4, 380)
        magnitude = offset + generator.normal(0, scatter, nights.size)
        for centre, depth, duration, _ in preregistered.KINDS:
            k = np.round((nights - centre) / preregistered.P26)
            magnitude += amplitude * preregistered.template(nights, centre + k * preregistered.P26, duration, depth)
        for i, (night, value) in enumerate(zip(nights, magnitude)):
            grade = "B" if i % 50 == 0 else "A"
            lines.append(f"{night + 7000:12.5f} {value:6.3f} {value:6.3f} {value:6.3f} {value:6.3f} {value:6.3f} 0.030 0.030 0.030 0.030 0.030  {grade} {i}")
        if other_block and dataset == 1:
            lines += [f"#dataset= {dataset} ; 2 F{dataset}", "#desig= 183416-6849.6", "-2000.00000  5.000  5.000  5.000  5.000  5.000 0.030 0.030 0.030 0.030 0.030  A 1".replace("-2000.00000", " 5000.00000")]
    path.write_text("\n".join(lines) + "\n")
    return path
