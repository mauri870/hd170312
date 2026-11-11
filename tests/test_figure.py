"""The figure script on synthetic sectors and a synthetic ASAS-3 file."""
import json

import numpy as np
import pytest
from synthetic_asas import synthetic_file
from synthetic_asassn import synthetic_file as asassn_file
from synthetic_neowise import EPOCH, target_file

from hd170312 import figure, neowise

SPANS = {13: (1653.9, 1682.4), 39: (2361.8, 2389.7), 66: (3097.7, 3121.6), 93: (3830.3, 3855.1), 94: (3857.0, 3881.1), 100: (4074.3, 4099.0), 101: (4101.4, 4125.3), 102: (4127.5, 4151.4), 103: (4153.3, 4177.4)}


def synthetic_sectors():
    """Nine sectors at the real spans, 20-minute cadence, flat with 300 ppm noise and the three complete eclipses and the end of sector 93 dipped."""
    generator = np.random.default_rng(0)
    sectors = {}
    for sector, (low, high) in SPANS.items():
        time = np.arange(low, high, 1 / 72)
        flux = 1 + generator.normal(0, 300e-6, time.size)
        for event in figure.EVENTS:
            if event["sector"] == sector:
                depth = 0.15 if event["kind"] == "deep" else 0.085  # the eclipse of sector 93 is in the SAP flux of the flagged cadences
                flux = flux * (1 - depth * np.clip((0.35 - abs(time - event["time"])) / 0.2, 0, 1))
        sectors[sector] = (time, flux)
    return sectors


def synthetic_sap():
    """The SAP flux of all nine sectors, so that no test reads the downloaded data: the sectors as in synthetic_sectors, the deep eclipse of sector 93 at its fitted time with its flagged cadences, ending 0.3 d after its minimum, and the shallow egress that begins sector 101."""
    generator = np.random.default_rng(1)
    out = {sector: {"time": time, "sap": flux * 7e4, "quality": np.zeros(time.size, dtype=int)} for sector, (time, flux) in synthetic_sectors().items() if sector not in (93, 101)}
    time = np.arange(3852.0, 3855.58, 1 / 720)
    flux = 1 + generator.normal(0, 300e-6, time.size)
    out[93] = {"time": time, "sap": flux * (1 - 0.155 * np.clip((0.33 - abs(time - 3855.2675)) / 0.06, 0, 1)) * 7e4, "quality": np.where(time > 3855.07, 4096, 0)}
    time = np.arange(4101.112, 4103.0, 1 / 720)
    flux = 1 + generator.normal(0, 300e-6, time.size)
    out[101] = {"time": time, "sap": flux * (1 - 0.085 * np.clip((0.37 - abs(time - 4100.8984)) / 0.08, 0, 1)) * 7e4, "quality": np.where(time < 4101.3, 4096, 0)}
    return out


def test_the_figure_is_drawn_and_panel_c_reports_the_nights_of_the_test(tmp_path, monkeypatch):
    monkeypatch.setattr(figure, "_SAP", synthetic_sap())
    out = tmp_path / "figure.png"
    nights, amplitude, error = figure.make(out, synthetic_file(tmp_path / "asas.txt"), neowise_path=target_file(tmp_path / "neowise.csv"), asassn_path=asassn_file(tmp_path / "asassn.csv"))
    assert out.stat().st_size > 50_000 and nights >= 6 and error > 0 and abs(amplitude - 1.0) < 0.5
    vector = out.with_suffix(".pdf")
    assert vector.exists() and vector.read_bytes()[:5] == b"%PDF-", "the figure is also written as a PDF beside the PNG, for the AAS"


def test_the_fourth_eclipse_is_a_deep_event_of_sector_93_in_the_sap_flux_and_the_egress_is_the_epoch_that_only_68_95_d_predicts():
    deep = [event for event in figure.EVENTS if event["sector"] == 93]
    assert len(deep) == 1 and deep[0]["kind"] == "deep" and abs(deep[0]["time"] - 3855.2675) < 1e-3 and deep[0]["source"] == "sap"
    assert abs(figure.EGRESS["time"] - (figure.SHALLOW_ANCHOR - figure.P26)) < 1e-3 and figure.EGRESS["sector"] == 101


def test_the_profile_of_the_flagged_eclipse_has_the_level_of_the_stretch_before_it_and_marks_its_flagged_cadences(monkeypatch):
    sap = synthetic_sap()
    sap[93]["sap"] = sap[93]["sap"] * (1 + 0.01 * (sap[93]["time"] - 3855.2675))  # a slope of the SAP flux, as in the real light curve
    monkeypatch.setattr(figure, "_SAP", sap)
    hours, relative, flagged = figure.profile(next(event for event in figure.EVENTS if event["sector"] == 93))
    before = (hours < -12) & (hours > -36)
    assert abs(np.median(relative[before]) - 1.0) < 1e-3 and flagged.sum() > 100 and flagged.size == hours.size
    assert not flagged[hours < -12].any(), "the cadences before the eclipse are not flagged"


def test_the_egress_of_panel_b_is_black_points_drawn_last_with_a_legend_handle_that_shows(monkeypatch):
    import matplotlib.pyplot as plt

    monkeypatch.setattr(figure, "_SAP", synthetic_sap())
    fig, ax = plt.subplots()
    figure.eclipses(ax)
    points = [line for line in ax.lines if line.get_marker() == "." and line.get_color() == "k"]
    assert len(points) == 1 and points[0].get_markersize() == 1.4 and points[0].get_alpha() == 0.45, "the earlier look: small black points"
    curves = [index for index, line in enumerate(ax.lines) if line.get_label().startswith("S")]
    assert ax.lines.index(points[0]) > max(curves), "drawn after the curves, so it stays on top"
    legend = ax.get_legend()
    entries = {text.get_text(): handle for text, handle in zip(legend.get_texts(), legend.legend_handles)}
    handle = entries["S101, possible egress (flagged)"]
    assert handle.get_linestyle() == "-" and handle.get_color() == "k" and handle.get_linewidth() == 1.2 and handle.get_marker() in ("None", None, ""), "a black line like the other entries, of their width"
    assert [text for text in entries if "egress" in text] == ["S101, possible egress (flagged)"], "one entry for the egress"
    plt.close(fig)


def test_the_four_eclipses_of_panel_b_are_drawn_from_the_sap_flux_and_the_figure_reads_no_pdcsap(monkeypatch):
    import matplotlib.pyplot as plt

    monkeypatch.setattr(figure, "_SAP", synthetic_sap())
    assert {event["source"] for event in figure.EVENTS} == {"sap"} and not hasattr(figure, "_SECTORS"), "PDCSAP is not read by the figure"
    event = next(event for event in figure.EVENTS if event["sector"] == 102)
    hours, relative, _ = figure.profile(event)
    assert abs(np.median(relative[(hours > 12) & (hours < 30)]) - np.median(relative[(hours < -12) & (hours > -30)])) < 1e-3, "no tilt across the eclipse"
    fig, ax = plt.subplots()
    figure.eclipses(ax)
    assert "S93 deep" in [text.get_text() for text in ax.get_legend().get_texts()] and "S93 deep (SAP)" not in [text.get_text() for text in ax.get_legend().get_texts()]
    plt.close(fig)


def test_the_lower_limit_of_panel_b_leaves_room_below_the_deepest_eclipse(monkeypatch):
    """The deep eclipse reaches about 0.845 in SAP flux; the axis must not cut it, so its lower limit is at least 0.01 below the deepest 1 h median of every eclipse."""
    import matplotlib.pyplot as plt

    monkeypatch.setattr(figure, "_SAP", synthetic_sap())
    fig, ax = plt.subplots()
    figure.eclipses(ax)
    deepest = min(np.nanmin(figure.hourly_medians(*figure.profile(event)[:2])[1]) for event in figure.EVENTS)
    low = ax.get_ylim()[0]
    plt.close(fig)
    assert low <= deepest - 0.01, f"lower limit {low} against the deepest median {deepest:.4f}"


def _panel_c(tmp_path, null):
    import matplotlib.pyplot as plt

    path = synthetic_file(tmp_path / "asas.txt")
    fig, ax = plt.subplots()
    errorbars = []
    original = ax.errorbar
    ax.errorbar = lambda *a, **k: (errorbars.append(k.get("yerr", a[2] if len(a) > 2 else None)), original(*a, **k))[1]
    numbers = figure.asas(ax, path, null=null)
    return fig, ax, path, errorbars, numbers


def test_the_error_bars_of_panel_c_are_the_errors_times_the_root_of_the_chi_square_per_degree_of_freedom_of_the_fit(tmp_path):
    from hd170312 import preregistered
    from hd170312.asas import nightly_means, read_rows, stacked

    path = synthetic_file(tmp_path / "asas.txt", scatter=0.030)  # noisier than the quoted errors, so that the scale is not 1
    arrays = stacked(nightly_means(read_rows(path)))
    (_, _, chi2, _), _ = preregistered.a_only(arrays)
    expected = np.sqrt(chi2 / (len(arrays[0]) - (len(set(arrays[3])) + 1)))
    assert figure.error_scale(arrays) == expected and expected > 1.02, "the scale is computed from the fit, not written in the code"
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    recorded = []
    original = ax.errorbar
    ax.errorbar = lambda *a, **k: (recorded.append(np.asarray(a[2])), original(*a, **k))[1]
    figure.asas(ax, path, null=np.linspace(0.0, 1.0, 200))
    plt.close(fig)
    errors = arrays[2]
    assert recorded and all(np.isin(np.round(bars / expected, 9), np.round(errors, 9)).all() for bars in recorded if bars.size), "every bar is a night's error times the scale"


def test_the_inset_of_panel_c_is_a_log_histogram_of_the_null_with_a_line_at_the_fitted_amplitude_and_its_p(tmp_path):
    import matplotlib.pyplot as plt

    null = np.random.default_rng(3).normal(0.0, 0.2, 2000)
    fig, ax, _, _, (_, amplitude, _) = _panel_c(tmp_path, null)
    inset = ax.child_axes[0]
    assert inset.get_yscale() == "log" and inset.get_xlabel() == "amplitude, shifted epochs"
    lines = [line for line in inset.lines if len(set(line.get_xdata())) == 1]
    assert len(lines) == 1 and lines[0].get_xdata()[0] == pytest.approx(amplitude)
    p = (1 + int((null >= amplitude).sum())) / (1 + null.size)
    assert any(f"p={p:.3f}" in text.get_text() and f"{amplitude:.2f}" in text.get_text() for text in inset.texts)
    plt.close(fig)


def test_the_null_is_read_from_the_stored_output_of_the_test_and_computed_only_when_there_is_none(tmp_path, monkeypatch):
    stored = tmp_path / "asas-preregistered.json"
    stored.write_text(json.dumps({"null_amplitudes": [0.1, 0.2, 0.3]}))
    assert list(figure.stored_null(stored)) == [0.1, 0.2, 0.3]
    assert figure.stored_null(tmp_path / "missing.json") is None
    stored.write_text(json.dumps({"p": 0.01}))
    assert figure.stored_null(stored) is None, "an output with no null amplitudes is not a null"


SPANS_ALL = {13: (1653.9, 1682.4), 39: (2361.8, 2389.7), 66: (3097.7, 3121.6), 93: (3830.3, 3855.6), 94: (3857.0, 3881.1), 100: (4074.3, 4100.5), 101: (4101.1, 4125.3), 102: (4127.5, 4151.4), 103: (4153.3, 4177.8)}


def sap_of_all_sectors():
    """The nine sectors with a cadence every 0.01 d, the last quarter of each flagged for scattered light."""
    out = {}
    for sector, (low, high) in SPANS_ALL.items():
        time = np.arange(low, high, 0.01)
        out[sector] = {"time": time, "sap": np.full(time.size, 7e4), "quality": np.where(time > low + 0.75 * (high - low), 4096, 0)}
    return out


def _phase(t):
    return ((np.asarray(t) - figure.SHALLOW_ANCHOR) / (figure.DELTA_T / 13)) % 1.0


def _panel_a(monkeypatch, tmp_path, arrays=None, sn_arrays=None):
    import matplotlib.pyplot as plt

    from hd170312.asas import nightly_means, read_rows, stacked

    monkeypatch.setattr(figure, "_SAP", {**synthetic_sap(), **sap_of_all_sectors(), 93: synthetic_sap()[93], 101: synthetic_sap()[101]})
    if arrays is None:
        arrays = stacked(nightly_means(read_rows(synthetic_file(tmp_path / "asas.txt"))))
    if sn_arrays is None:
        sn_arrays = _nights_of_the_test()
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    figure.phase_diagram(ax, arrays, neowise.read_exposures(target_file(tmp_path / "neowise.csv")), sn_arrays)
    return fig, ax


def test_the_four_lines_of_panel_a_are_at_the_phases_of_the_anchors_over_the_137_90_d_orbit_and_half_an_orbit_later():
    deep = ((figure.DEEP_ANCHOR - figure.SHALLOW_ANCHOR) / (figure.DELTA_T / 13)) % 1.0
    assert [(round(x, 6), kind, only) for x, kind, only in figure.LINES] == [(0.0, "shallow", False), (round(deep, 6), "deep", False), (0.5, "shallow", True), (round((deep + 0.5) % 1.0, 6), "deep", True)]
    assert abs(deep - 0.71876) < 1e-4


def test_every_observed_eclipse_lies_on_a_solid_line_and_the_egress_and_the_2027_eclipse_on_the_dashed_ones():
    solid = {kind: x for x, kind, only in figure.LINES if not only}
    dashed = {kind: x for x, kind, only in figure.LINES if only}
    for event in figure.EVENTS:
        assert abs(figure.marker_phase(event["time"]) - solid[event["kind"]]) < 0.01, event
    assert figure.marker_phase(figure.EVENTS[-1]["time"]) < 0.01 and figure.marker_phase(figure.EVENTS[-1]["time"]) > -0.01, "the eclipse of sector 103 is at phase 0, wrapped, not 1"
    assert abs(figure.marker_phase(figure.EGRESS["time"]) - dashed["shallow"]) < 0.01
    assert abs(figure.marker_phase(figure.decisive_eclipse()) - dashed["deep"]) < 0.01 and abs(figure.decisive_eclipse() - 4475.815) < 0.01


def test_the_bar_of_a_sector_covers_the_phases_of_its_cadences_in_dark_grey_for_good_ones_and_light_grey_for_flagged_ones(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = _panel_a(monkeypatch, tmp_path)
    sap = figure._SAP
    for sector in (13, 66, 101):
        y = figure.row_y(sector)
        bars = [patch for patch in ax.patches if abs(patch.get_y() + patch.get_height() / 2 - y) < 1e-9]
        for flagged, colour in ((False, "0.55"), (True, "0.8")):
            mask = (sap[sector]["quality"] == 4096) == flagged
            phases = _phase(sap[sector]["time"][mask])
            mine = [(patch.get_x(), patch.get_x() + patch.get_width()) for patch in bars if patch.get_facecolor()[:3] == matplotlib_color(colour)]
            if not phases.size:
                continue
            assert all(any(a - 1e-9 <= x <= b + 1e-9 for a, b in mine) for x in phases), f"S{sector} {'flagged' if flagged else 'good'} cadences are inside its bars"
            assert all(any(abs(x - (a + b) / 2) <= (b - a) / 2 + figure.PHASE_STEP for x in phases) for a, b in mine), "no bar covers phases that have no cadence"
    plt.close(fig)


def matplotlib_color(grey):
    from matplotlib.colors import to_rgb

    return to_rgb(grey)


def test_the_lines_of_panel_a_stop_below_their_labels(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = _panel_a(monkeypatch, tmp_path)
    labels = [text for text in ax.texts if "both periods" in text.get_text() or "68.95 d only" in text.get_text()]
    lines = [line for line in ax.lines if len(set(line.get_xdata())) == 1 and len(line.get_ydata()) == 2]
    assert len(labels) == 4 and len(lines) == 4
    assert all(max(line.get_ydata()) < min(text.get_position()[1] for text in labels) for line in lines), "a line ends under the labels"
    plt.close(fig)


def test_every_symbol_drawn_in_panel_a_has_an_entry_in_its_legend(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = _panel_a(monkeypatch, tmp_path)

    from matplotlib.colors import to_rgb

    def signature(handle):
        if hasattr(handle, "get_marker") and handle.get_marker() == "|":
            return ("tick", tuple(round(v, 3) for v in to_rgb(handle.get_color())))
        if hasattr(handle, "get_marker"):
            return ("marker", handle.get_marker(), handle.get_markerfacecolor() in ("white", "none"), str(handle.get_markeredgecolor() if handle.get_markerfacecolor() in ("white", "none") else handle.get_markerfacecolor()))
        return ("patch", handle.get_hatch() or "", tuple(round(v, 3) for v in handle.get_facecolor()[:3]) if not handle.get_hatch() else "hatched")

    drawn = {signature(line) for line in ax.lines if line.get_marker() not in ("None", None, "")}
    drawn |= {signature(patch) for patch in ax.patches}
    drawn |= {("tick", tuple(round(v, 3) for v in collection.get_colors()[0][:3])) for collection in ax.collections}
    legend = {signature(handle) for handle in ax.get_legend().legend_handles}
    assert drawn <= legend, drawn - legend
    plt.close(fig)


def _nights_of_the_test():
    """Thirty nights: ten at odd epochs of 68.95 d (five of each kind, inside half a duration, off every even epoch) and twenty far from every eclipse."""
    from hd170312 import preregistered

    inside = [centre + k * preregistered.P26 + 0.05 for centre, _, _, _ in preregistered.KINDS for k in (-33, -31, -29, -27, -25)]
    far = [2500.0 + 7.3 * i for i in range(20)]
    T = np.sort(np.array(inside + far))
    return T, np.full(T.size, 8.7), np.full(T.size, 0.01), np.ones(T.size, dtype=int), np.ones(T.size)


def test_the_orange_ticks_of_the_asas_row_are_exactly_the_nights_of_the_test_and_the_grey_ticks_are_the_others(monkeypatch, tmp_path):
    from hd170312 import preregistered

    arrays = _nights_of_the_test()
    tested = preregistered.only_window(arrays[0])
    assert tested.sum() == 10, "the nights of the test: ten, as in the real data"
    fig, ax = _panel_a(monkeypatch, tmp_path, arrays)
    ticks = {tuple(round(v, 3) for v in collection.get_colors()[0][:3]): sorted(segment[0][0] for segment in collection.get_segments()) for collection in ax.collections}
    from matplotlib.colors import to_rgb

    orange, grey = ticks[tuple(round(v, 3) for v in to_rgb(figure.ORANGE))], ticks[tuple(round(v, 3) for v in to_rgb("0.55"))]
    assert np.allclose(orange, sorted(figure.phase(arrays[0][tested])), atol=1e-12) and len(orange) == 10
    assert len(grey) == arrays[0].size - 10 and np.allclose(grey, sorted(figure.phase(arrays[0][~tested])), atol=1e-12)
    # every orange tick lies within half an eclipse duration of an orange line
    half = max(duration for _, _, duration, _ in preregistered.KINDS) / 2 / figure.P13
    lines = [x for x, _, only in figure.LINES if only]
    assert all(min(abs(x - line) for line in lines) <= half + 1e-9 for x in orange), "each tick is inside an eclipse that only 68.95 d predicts"
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_the_neowise_row_is_between_sector_13_and_the_asas_row_and_the_labels_of_the_lines_are_above_them(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = _panel_a(monkeypatch, tmp_path)
    assert next(text for text in ax.texts if text.get_text() == "NEOWISE (2014–24)").get_position()[1] == figure.row_y(13) + 1
    label = next(text for text in ax.texts if text.get_text() == "ASAS-3 (2000–09)")
    assert label.get_position()[1] == figure.row_y(13) + 2
    assert all(text.get_position()[1] > label.get_position()[1] for text in ax.texts if text.get_text().startswith(("both", "68.95")))
    plt.close(fig)


def test_each_visit_of_neowise_is_a_grey_bar_at_its_phase_and_one_of_them_lies_on_the_dashed_deep_line(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb

    fig, ax = _panel_a(monkeypatch, tmp_path)
    rows = neowise.read_exposures(target_file(tmp_path / "neowise.csv"))
    y = figure.row_y(13) + 1
    bars = [patch for patch in ax.patches if abs(patch.get_y() + patch.get_height() / 2 - y) < 1e-9]
    grey = [patch for patch in bars if patch.get_facecolor()[:3] == to_rgb(figure.NEOWISE_GREY)]
    assert len(bars) == 28 and len(grey) == 28, "28 visits, all grey; the triangle marks the one with the eclipse"
    deep_only = next(x for x, kind, only in figure.LINES if only and kind == "deep")
    assert sum(patch.get_x() <= deep_only <= patch.get_x() + patch.get_width() for patch in grey) == 1, "one visit lies on the dashed deep line"
    assert abs(figure.phase(EPOCH) - deep_only) < 0.01
    visit_phases = [figure.phase(rows["btjd"][rows["visit"] == v]) for v in range(28)]
    for patch in grey:
        assert any(patch.get_x() - 1e-9 <= p <= patch.get_x() + patch.get_width() + 1e-9 for phases in visit_phases for p in phases)
    labels = [handle.get_label() for handle in ax.get_legend().legend_handles]
    assert "NEOWISE visit" in labels and "NEOWISE visit with the eclipse only 68.95 d predicts" in labels
    plt.close(fig)


def test_the_neowise_visit_with_the_eclipse_that_only_68_95_d_predicts_carries_a_triangle_so_it_is_not_read_as_an_asas_3_tick(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = _panel_a(monkeypatch, tmp_path)
    y = figure.row_y(13) + 1
    triangles = [line for line in ax.lines if line.get_marker() == "v" and abs(line.get_ydata()[0] - y) < 1e-9]
    assert len(triangles) == 1, "one triangle on the NEOWISE row"
    deep_only = next(x for x, kind, only in figure.LINES if only and kind == "deep")
    assert abs(triangles[0].get_xdata()[0] - deep_only) < 0.01, "the triangle sits on the dashed deep line"
    handle = next(handle for handle in ax.get_legend().legend_handles if handle.get_label() == "NEOWISE visit with the eclipse only 68.95 d predicts")
    assert handle.get_marker() == "v", "the legend shows the triangle"
    plt.close(fig)


def test_the_asassn_row_is_above_the_asas_3_row_with_every_night_a_light_grey_tick_and_the_nights_of_the_windows_hollow_orange_rectangles(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb

    from hd170312 import preregistered

    sn = _nights_of_the_test()
    tested = preregistered.only_window(sn[0])
    fig, ax = _panel_a(monkeypatch, tmp_path, sn_arrays=sn)
    label = next(text for text in ax.texts if text.get_text() == "ASAS-SN (2014–26)")
    assert label.get_position()[1] == figure.row_y(13) + 3 and next(text for text in ax.texts if text.get_text() == "ASAS-3 (2000–09)").get_position()[1] == figure.row_y(13) + 2
    grey = [c for c in ax.collections if tuple(round(v, 3) for v in c.get_colors()[0][:3]) == tuple(round(v, 3) for v in to_rgb(figure.SN_GREY))]
    assert len(grey) == 1 and np.allclose(sorted(segment[0][0] for segment in grey[0].get_segments()), sorted(figure.phase(sn[0])), atol=1e-12), "every night is a light grey tick at its phase"
    hollow = [patch for patch in ax.patches if abs(patch.get_y() + patch.get_height() / 2 - label.get_position()[1]) < 1e-9]
    assert len(hollow) == int(tested.sum()) == 10
    assert all(not patch.get_fill() and patch.get_edgecolor()[:3] == to_rgb(figure.ORANGE) for patch in hollow), "outline only, orange"
    dashed = {kind: x for x, kind, only in figure.LINES if only}
    for centre_time, _, duration, kind in preregistered.KINDS:
        boxes = [patch for patch in hollow if abs(patch.get_width() - duration / figure.P13) < 1e-9]
        assert len(boxes) == 5 and all(abs(patch.get_x() + patch.get_width() / 2 - dashed[kind]) < 1e-6 for patch in boxes), f"five {kind} windows, half a duration each side of the epoch on its dashed line"
    labels = [handle.get_label() for handle in ax.get_legend().legend_handles]
    assert "ASAS-SN night" in labels and "ASAS-SN night in an eclipse only 68.95 d predicts" in labels
    assert all(patch.get_linewidth() < 1.3 for patch in hollow), "thinner than the filled ASAS-3 ticks (1.3)"
    plt.close(fig)


def test_the_figure_is_drawn_with_the_asassn_file(tmp_path, monkeypatch):
    monkeypatch.setattr(figure, "_SAP", synthetic_sap())
    out = tmp_path / "figure.png"
    figure.make(out, synthetic_file(tmp_path / "asas.txt"), neowise_path=target_file(tmp_path / "neowise.csv"), asassn_path=asassn_file(tmp_path / "asassn.csv"))
    assert out.stat().st_size > 50_000


def test_the_row_labels_of_the_sectors_name_the_instrument(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = _panel_a(monkeypatch, tmp_path)
    texts = {text.get_text() for text in ax.texts}
    assert {"TESS S13 (2019)", "TESS S103 (2026)", "TESS S114 (2027, planned)"} <= texts and not any(t.startswith("S1") for t in texts)
    plt.close(fig)


def test_a_bar_runs_from_the_phase_of_its_first_cadence_to_that_of_its_last_and_not_out_to_the_bins_of_the_phase_step(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = _panel_a(monkeypatch, tmp_path)
    sap = figure._SAP
    y = figure.row_y(101)
    bars = [patch for patch in ax.patches if abs(patch.get_y() + patch.get_height() / 2 - y) < 1e-9]
    first = float(_phase(sap[101]["time"][0]))
    assert abs(first - 0.5015) < 2e-4, "the first cadence of sector 101 is at BTJD 4101.112, phase 0.5015, after the egress marker at 0.5"
    assert abs(min(patch.get_x() for patch in bars) - first) < 1e-4, "the bar of sector 101 starts at its first cadence, not under the egress marker at 0.500"
    last = float(_phase(sap[101]["time"][-1]))
    assert abs(max(patch.get_x() + patch.get_width() for patch in bars) - last) < 1e-4
    plt.close(fig)


def test_segments_split_where_the_cadences_leave_a_gap_longer_than_the_step_and_keep_a_single_cadence_as_a_point():
    times = figure.SHALLOW_ANCHOR + figure.P13 * np.array([0.1000, 0.1005, 0.1010, 0.3000, 0.3003, 0.8])
    found = figure.segments(times)
    assert [(round(a, 4), round(b, 4)) for a, b in found] == [(0.1, 0.101), (0.3, 0.3003), (0.8, 0.8)]
    assert figure.segments(np.array([])) == []


def test_the_legend_calls_the_egress_possible_and_flagged_and_the_2027_eclipse_predicted(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = _panel_a(monkeypatch, tmp_path)
    labels = [handle.get_label() for handle in ax.get_legend().legend_handles]
    assert "possible egress, flagged data (68.95 d only)" in labels and "predicted eclipse, 2027 March 11" in labels
    assert not any(label in ("egress seen (68.95 d only)", "2027 March 11") for label in labels)
    plt.close(fig)
