import pathlib

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca
from development.time_series_support import ResultTimeSeries, plot_results

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library_color.json"
LIBRARY = ca.gc_lc.GCLibrary.from_JSON(lib_path)
INTERNAL_STANDARD = LIBRARY.find_by_label("TCB")
picking_lib_fid = ca.a.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_fid")
picking_lib_ms = ca.a.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_ms")
PLOTTING_GROUPS = ["dicarboxylic acid", "hydroxy acids", "carboxylic acid", "alcohol", "methyl_ketone", "ketone",
                   "peroxide", "alkane"]


def process_fid(fid: ca.gc_lc.GCSignal2D):
    figs = []

    # processing
    fid.processor.add(
        ca.p.edit.ReplaceSpans(value=0,
                               x_spans=(
                                   (1.3, 3.7),  # solvent
                                   (8.7, 10),  # decane
                                   (44.15, 44.85)  # PPh3
                               ), invert=True),
        ca.p.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4)
    )

    # analysis
    results = ResultTimeSeries("decane_fid", INTERNAL_STANDARD, internal_standard_mmol=0.0109)
    for sig in fid.signal_iter():
        peak_locations = ca.a.peak_picking.find_peaks_scipy(sig, scipy_kwargs={"height": 8000, "width": 0.1})
        peaks = ca.a.integration.rolling_ball(peak_locations, n=5, min_height=0.002, n_points_with_pos_slope=2)
        compounds = ca.a.ms_analysis.search_by_retention_time(picking_lib_fid, peaks)

        fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
        ca.plotting.signal(sig, fig=fig)
        ca.plotting.peaks(compounds, fig=fig)
        figs.append(fig)

        results.add_result(compounds, sig.extract_value)

    # plotting results
    # fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    # plot_results(fid_timeseries, PLOTTING_GROUPS, fig=fig)
    # # fig.add_scatter(x=[0, 360], y=[0.0658, 0.0658], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    # fig.layout.xaxis.title = "<b>time (min)<br>"
    # fig.layout.yaxis.title = "<b>mmol<br>"
    # figs.append(fig)
    #
    # fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    # plot_results(fid_timeseries, PLOTTING_GROUPS, fig=fig, carbon=True)
    # fig.add_scatter(x=[0, fid.y[-1]], y=[0.0658*10, 0.0658*10], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    # fig.layout.xaxis.title = "<b>time (min)<br>"
    # fig.layout.yaxis.title = "<b>mmol of carbon<br>"
    # figs.append(fig)

    figs_html = ca.plotting.plotly_utils.merge_figures(figs, auto_open=True)
    return figs_html


def process_ms(ms: ca.gc_lc.GCMSSignal2D):
    figs = []

    # processing
    ms.processor.add(
        ca.p.edit.ReplaceSpans(value=0,
                       x_spans=(
                           (8.9, 10),  # decane
                           (44.1, 45.1)  # PPh3
                       ), invert=True),
        ca.p.baseline.SectionMinMax(sections=100, window=50, number_of_deviations=4,
                                    smoother=ca.p.smoothing.Gaussian())
    )

    # analysis
    results = ResultTimeSeries("decane_ms", INTERNAL_STANDARD, internal_standard_mmol=0.0109)
    for sig in ms.signal_iter():
        peak_locations = ca.a.peak_picking.find_peaks_scipy(sig, scipy_kwargs={"prominence": 10000})
        peaks = ca.a.integration.rolling_ball(peak_locations, n=5, min_height=0.002, n_points_with_pos_slope=2)
        compounds = ca.a.ms_analysis.search_by_retention_time(picking_lib_ms, peaks)

        fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
        ca.plotting.signal(sig, fig=fig)
        ca.plotting.peaks(compounds, fig=fig)
        figs.append(fig)

        results.add_result(compounds, sig.extract_value)

    # fig = go.Figure(layout=ca.p.plotly_utils.layout())
    # plot_results(results, PLOTTING_GROUPS, fig=fig)
    # # fig.add_scatter(x=[0, 360], y=[0.0658, 0.0658], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    # fig.layout.xaxis.title = "<b>time (min)<br>"
    # fig.layout.yaxis.title = "<b>mmol<br>"
    # figs.append(fig)
    #
    # fig = go.Figure(layout=ca.p.plotly_utils.layout())
    # plot_results(results, PLOTTING_GROUPS, fig=fig, carbon=True)
    # fig.add_scatter(x=[0, ms.y[-1]], y=[0.0658 * 10, 0.0658 * 10], mode="lines",
    #                 line={"color": "black", "dash": "dash"}, name="decane_init")
    # fig.layout.xaxis.title = "<b>time (min)<br>"
    # fig.layout.yaxis.title = "<b>mmol of carbon<br>"
    # figs.append(fig)

    data = results.to_csv_str()
    print(data)
    figs_html = ca.plotting.plotly_utils.merge_figures(figs, auto_open=True)
    return figs_html


def main():
    data_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\10\10_10")
    fid_timeseries = ca.gc_lc.GCSignal2D.from_npz(data_path / "FID_timeseries.npz")
    ms_timeseries = ca.gc_lc.GCMSSignal2D.from_npz(data_path / "MS_timeseries.npz")
    ms_figs = process_ms(ms_timeseries)
    fid_figs = process_fid(fid_timeseries)


if __name__ == "__main__":
    main()
