"""A lower bound on the radius of the body B that eclipses the star A, from the depth of an eclipse: B in front of A blocks the share (R_B / R_A)^2 of the light of the disk it covers.

Without limb darkening the depth is (R_B / R_A)^2. With it, an eclipse across the centre of the disk, which is the smallest body for a depth, hides light that is I_centre / mean(I) times brighter than the mean of the disk, so
R_B >= R_A sqrt(depth / (I_centre / mean(I))); a graze or the light of B only raise the bound. For the quadratic law I(mu) / I(1) = 1 - a (1 - mu) - b (1 - mu)^2 the mean of the disk is I(1) (1 - a / 3 - b / 6)."""
import numpy as np

# TESS-band quadratic coefficients (aLSM, bLSM) of Claret (2017, A&A 600, A30), ATLAS models, [M/H] = 0, microturbulence 2 km/s, log g 3.5: table25 of https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/600/A30
# (VizieR J/A+A/600/A30/table25, read on 2026-11-11 at the mirror vizier.nao.ac.jp with Teff = 6500 to 6750, logg = 3.5, Z = 0, xi = 2: aLSM 0.2028 and 0.1843, bLSM 0.3226 and 0.3280).
CLARET = {6500.0: (0.2028, 0.3226), 6750.0: (0.1843, 0.3280)}
TEFF = 6645.0  # K, TESS Input Catalog v8.2 (the effective temperature of the SPOC data validation report of the star, whose log g is 3.51)
# The radii of A in solar radii: the TIC radius (Stassun et al. 2019, AJ 158, 138) and the radii that fits of two spherical stars need at 68.95 d and 137.90 d, which the Note gives (24 and 34 percent below the TIC radius);
# the repository does not compute the fits of two stars.
RADII_OF_A = {"TIC": 3.46, "two stars, 68.95 d": 2.64, "two stars, 137.90 d": 2.27}


def limb_darkening(teff):
    """(a, b) at a temperature between the two rows of CLARET, interpolated linearly in Teff; ValueError outside them."""
    low, high = sorted(CLARET)
    if not low <= teff <= high:
        raise ValueError(f"Teff {teff} is outside {low:.0f} to {high:.0f} K, the two rows of the table of Claret (2017) that this module holds")
    fraction = (teff - low) / (high - low)
    return tuple(float(CLARET[low][index] + fraction * (CLARET[high][index] - CLARET[low][index])) for index in (0, 1))


def central_over_mean(a, b):
    """I_centre / mean(I) of a disk with the quadratic law of coefficients a and b."""
    return 1.0 / (1.0 - a / 3.0 - b / 6.0)


def lower_bound(depth, radius_of_a, coefficients=None):
    """The radius of B in solar radii: sqrt(depth) times the radius of A, over the square root of I_centre / mean(I) when the coefficients (a, b) are given."""
    if depth <= 0:
        raise ValueError(f"a depth of {depth} gives no bound")
    brightness = 1.0 if coefficients is None else central_over_mean(*coefficients)
    return float(np.sqrt(depth / brightness) * radius_of_a)


def table(events):
    """{limb_darkening, cases, conservative}: for each kind of eclipse (its depth is that of the shallower eclipse of the kind) and each radius of A, the bound without and with limb darkening; `conservative` is the smallest of all.

    events are the fits of hd170312.events. The bound holds whichever eclipse is B in front of A: the shallow eclipse gives the smaller one."""
    coefficients = limb_darkening(TEFF)
    depths = {kind: min(event["depth"] for event in events.values() if event["kind"] == kind) for kind in ("shallow", "deep")}
    cases = [{"kind": kind, "depth": depths[kind], "radius_of_a_from": source, "radius_of_a": radius_of_a, "limb_darkening": limb, "radius_of_b": lower_bound(depths[kind], radius_of_a, coefficients if limb else None)}
             for kind in depths for source, radius_of_a in RADII_OF_A.items() for limb in (False, True)]
    return {"limb_darkening": {"teff": TEFF, "a": coefficients[0], "b": coefficients[1], "central_over_mean": central_over_mean(*coefficients)}, "cases": cases, "conservative": min(cases, key=lambda case: case["radius_of_b"])}
