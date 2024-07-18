import pathlib

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca
from development.time_series_support import ResultTimeSeries, plot_results

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library_color.json"
LIBRARY = ca.gc_lc.GCLibrary.from_JSON(lib_path)
INTERNAL_STANDARD = LIBRARY.find_by_label("TCB")
picking_lib_fid = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_fid")
picking_lib_ms = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_ms")
PLOTTING_GROUPS = ["dicarboxylic acid", "hydroxy acids", "carboxylic acid", "alcohol", "methyl_ketone", "ketone", "peroxide", "alkane"]


def process_fid(fid: ca.gc_lc.GCSignal2D):
    # processing
    fid.processor.add(
        ca.processing.edit.ReplaceSpans(value=0,
                                        x_spans=(
                                            (1.3, 3.7),  # solvent
                                            (8.7, 10),  # decane
                                            (44.15, 44.85)  # PPh3
                                        ), invert=True),
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
    )

    # # analysis
    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(fid, scipy_kwargs={"height": 8000, "width": 0.1})
    fid_peaks = ca.analysis.integration.rolling_ball(peak_locations, n=5, min_height=0.002,
                                                     n_points_with_pos_slope=2)
    # fid_compounds = ca.analysis.ms_analysis.search_by_retention_time(picking_lib_fid, fid_peaks)
    #
    # fid_timeseries = ResultTimeSeries("decane_fid", INTERNAL_STANDARD, internal_standard_mmol=0.0109)
    # for i in range(len(fid.y)):
    #     fid_timeseries.add_result(fid_compounds[i], fid.y[i])

    # plotting
    figs = []
    fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    figs_ = ca.plotting.signal2D_slices_separate(fid, fig=fig)
    ca.plotting.signal2D_slices_peaks(fid_peaks, fig=figs_)
    ca.plotting.plotly_utils.merge_figures(figs_, auto_open=True)
    exit()
    figs.append(fig)

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

    figs_html = ca.plotting.PlotlyConfig.merge_figures(figs)
    return figs_html

# def process_ms(ms: ca.gc_lc.GCMSSignal2D):
#     ms.processor.add(
#         ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
#     )
#     peak_locations = ca.analysis.peak_picking.find_peaks_scipy(ms,
#                                                                scipy_kwargs={"prominence": 10000}
#                                                                )
#     ms_peaks = ca.analysis.integration.rolling_ball(peak_locations, n=5, min_height=0.002,
#                                                     n_points_with_pos_slope=2)
#     ms_compounds = ca.analysis.ms_analysis.search_by_retention_time(picking_lib_ms, ms_peaks)
#
#     # plotting peak results
#     ms_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
#     ca.plotting.signal(ms, fig=ms_fig)
#     ca.plotting.peaks(ms_compounds, fig=ms_fig)
#     ms_fig.layout.title = "MS"
#
#     ms_timeseries = ResultTimeSeries("decane_ms", INTERNAL_STANDARD, internal_standard_mmol=0.0109)
#     for i in range(len(times)):
#         ms_timeseries.add_result(ms_compounds[i], times[i])
#
#         fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
#     plot_results(ms_timeseries, PLOTTING_GROUPS, fig=fig)
#     # fig.add_scatter(x=[0, 360], y=[0.0658, 0.0658], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
#     fig.layout.xaxis.title = "<b>time (min)<br>"
#     fig.layout.yaxis.title = "<b>mmol<br>"
#     ms_figs.append(fig)
#
#         fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
#     plot_results(ms_timeseries, PLOTTING_GROUPS, fig=fig, carbon=True)
#     fig.add_scatter(x=[0, times[-1]], y=[0.0658*10, 0.0658*10], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
#     fig.layout.xaxis.title = "<b>time (min)<br>"
#     fig.layout.yaxis.title = "<b>mmol of carbon<br>"
#     ms_figs.append(fig)
#
#     ca.plotting.PlotlyConfig.merge_figures(ms_figs, filename=figure_folder / "ms")
#     data = ms_timeseries.to_csv_str()
#     print(data)

def main():
    data_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\10\10_10")
    fid_timeseries = ca.gc_lc.GCSignal2D.from_npz(data_path / "FID_timeseries.npz")
    ms_timeseries = ca.gc_lc.GCMSSignal2D.from_npz(data_path / "MS_timeseries.npz")
    # process_ms(ms_timeseries)
    fid_figs = process_fid(fid_timeseries)

    # saving
    with open(data_path / "FID_timeseries.html", "w", encoding="UTF-8") as f:
        f.write(fid_figs)

if __name__ == "__main__":
    main()
