"""Synthetic IRSA TAP results for the offline tests: the 10 arcsec cone of the star (the columns of data.NEOWISE_COLUMNS) and the 600 arcsec field of the visit with the eclipse that only 68.95 d predicts. The real files are fetched, never committed."""
import numpy as np

from hd170312 import neowise, preregistered

TARGET_HEADER = "ra,dec,mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro,cc_flags,ph_qual,qual_frame"
FIELD_HEADER = "ra,dec,mjd,w1mpro,w2mpro,cc_flags,ph_qual,qual_frame,allwise_cntr,w2mpro_allwise"
EPOCH = -74.870  # the deep eclipse of the visit that only 68.95 d predicts, BTJD


def mjd_of(btjd):
    """The mjd whose BTJD is `btjd`."""
    mjd = btjd + 56999.5
    return mjd + (btjd - float(neowise.btjd_of(np.array([mjd]))[0]))


def visit_times():
    """Twenty-eight visits of 13 exposures 0.1 d apart in time order: the third visit (index 2) centred on EPOCH, the others at least 3 d from every epoch of either kind."""
    starts = [-260.0, -170.0]
    candidate = 100.0
    while len(starts) < 27:
        if neowise.distance_to_epoch(candidate + 0.1 * np.arange(13)).min() > 3.0:
            starts.append(candidate)
        candidate += 37.1
    visits = [start + 0.1 * np.arange(13) for start in starts]
    visits.insert(2, EPOCH + 0.1 * (np.arange(13) - 6))
    return visits


def target_file(path, dip=0.1536, seed=3):
    """Write the cone of the star: 28 visits, W2 flat at 7.57 with the deep template of flux depth `dip` in the visit of EPOCH, 0.02 mag of noise, and one row for each cut: qual_frame 0, a flagged cc_flags, ph_qual AX, a source 6 arcsec away."""
    generator = np.random.default_rng(seed)
    lines = [TARGET_HEADER]
    for btjd in visit_times():
        for t in btjd:
            w2 = 7.57 + float(preregistered.template(np.array([t]), EPOCH, 0.653, dip)[0]) + generator.normal(0, 0.02)
            lines.append(f"278.5782800,-68.8325300,{mjd_of(t):.8f},{7.45 + generator.normal(0, 0.06):.3f},0.018,{w2:.3f},0.013,0000,AA,10")
    time = mjd_of(btjd[0])
    lines += [f"278.5782800,-68.8325300,{time:.8f},7.45,0.018,7.57,0.013,0000,AA,0", f"278.5782800,-68.8325300,{time:.8f},7.45,0.018,7.57,0.013,dH00,AA,10", f"278.5782800,-68.8325300,{time:.8f},7.45,0.018,7.57,0.013,0000,AX,10",
              f"278.5782800,{-68.8325300 + 6 / 3600:.7f},{time:.8f},7.45,0.018,7.57,0.013,0000,AA,10"]
    path.write_text("\n".join(lines) + "\n")
    return path


def field_file(path, shifts, seed=5):
    """Write the field of the 13 frames of the visit of EPOCH: the star itself and one comparison star for each value of `shifts` (mag fainter in the window frames than in the others), AllWISE W2 8 to 10, 300 arcsec from the star; and one source with no AllWISE counterpart, one fainter than W2 11 and one 5 arcsec from the star, which are not comparison stars."""
    generator = np.random.default_rng(seed)
    times = visit_times()[2]
    window = preregistered.only_window(times)
    lines = [FIELD_HEADER]
    stars = [("1", 278.5782800, -68.8325300, 7.57, 0.0, 7.58)]
    stars += [(str(10 + i), 278.5782800 + 0.02 * (i + 1), -68.8325300, 8.0 + 0.2 * i, shift, 8.0 + 0.2 * i) for i, shift in enumerate(shifts)]
    stars += [("", 278.60, -68.83, 9.0, 0.0, 9.0), ("99", 278.62, -68.84, 12.5, 0.0, 12.5), ("98", 278.5782800, -68.8325300 + 5 / 3600, 9.0, 0.0, 9.0)]
    for cntr, ra, dec, w2, shift, catalogue in stars:
        for t, inside in zip(times, window):
            lines.append(f"{ra:.7f},{dec:.7f},{mjd_of(t):.8f},{w2 + 0.1:.3f},{w2 + shift * inside + generator.normal(0, 0.003):.3f},0000,AA,10,{cntr},{catalogue:.3f}")
    path.write_text("\n".join(lines) + "\n")
    return path
