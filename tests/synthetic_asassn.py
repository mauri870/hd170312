"""A synthetic ASAS-SN Sky Patrol v1 file in the layout of the saturated-star photometry download (a header line, a comment line of camera offsets, one row for each measurement), for the offline tests: the real file is tracked in data/asassn."""
import numpy as np

from hd170312 import preregistered

HEADER = "HJD,UT Date,Camera,FWHM,Limit,mag,mag_err,flux(mJy),flux_err,Filter"


def synthetic_file(path, odd=1.0, even=2.0, seed=3, scatter=0.010, per_night=2):
    """Write the file: two cameras (g:bi and V:be) with `per_night` points a night on every night of 2400 d, the eclipses of the Note injected at `odd` times the TESS depth on the epochs only 68.95 d predicts and at `even` times on the epochs both periods predict; a row with no magnitude is added. Returns the path."""
    generator = np.random.default_rng(seed)
    lines = [HEADER, "### magnitude offsets  bi -0.000  be -0.264 "]
    for camera, band, offset in (("bi", "g", 8.70), ("be", "V", 8.90)):
        nights = np.arange(-3000, -600, 1.0) + 0.3
        for night in nights:
            for point in range(per_night):
                hjd = night + 0.01 * point
                magnitude = offset + generator.normal(0, scatter)
                for centre, depth, duration, _ in preregistered.KINDS:
                    k = np.round((hjd - centre) / preregistered.P26)
                    amplitude = odd if k % 2 != 0 else even
                    magnitude += amplitude * preregistered.template(np.array([hjd]), centre + k * preregistered.P26, duration, depth)[0]
                lines.append(f"{hjd + 2457000:.6f},2014-1-1.0,{camera},2.0,NA,{magnitude:.3f},0.030,NA,NA,{band}")
    lines.append("2456800.500000,2014-5-17.0,bi,2.0,NA,NA,NA,NA,NA,g")
    path.write_text("\n".join(lines) + "\n")
    return path
