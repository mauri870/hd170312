"""The three complete eclipses of HD 170312 in the PDCSAP flux of sectors 39, 102 and 103."""
from hd170312 import data
from hd170312.eclipses import fit_event_time

EVENTS = [("S39", 39, 2377.16, "shallow"), ("S102", 102, 4131.07, "deep"), ("S103", 103, 4169.86, "shallow")]  # name, sector, a time within 0.1 d of the eclipse, kind


def fit_events(directory=data.DATA, loader=data.load_sector):
    """{name: the fit of fit_event_time plus its sector and kind} for the three complete eclipses."""
    found = {}
    for name, sector, guess, kind in EVENTS:
        columns = loader(sector, directory) if loader is data.load_sector else loader(sector)
        fit = fit_event_time(columns["time"], columns["pdcsap"], guess)
        found[name] = {**fit, "sector": sector, "kind": kind}
    return found
