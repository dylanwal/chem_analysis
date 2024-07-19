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


def process_single(data_path: pathlib.Path, label: str):
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)

    # fid
    fid.processor.add(
        ca.p.edit.ReplaceSpans(value=0,
                               x_spans=(
                                   (1.3, 3.7),  # solvent
                                   (8.7, 10),  # decane
                                   (44.15, 44.85)  # PPh3
                               ), invert=True),
        ca.p.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4)
    )
    peak_locations = ca.a.peak_picking.find_peaks_scipy(fid, scipy_kwargs={"height": 8000, "width": 0.1})
    peaks = ca.a.integration.rolling_ball(peak_locations, n=5, min_height=0.002, n_points_with_pos_slope=2)
    fid_compounds = ca.a.ms_analysis.search_by_retention_time(picking_lib_fid, peaks)
    
    fid_fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    ca.plotting.signal(fid, fig=fid_fig)
    ca.plotting.peaks(fid_compounds, fig=fid_fig)
    fid_fig.layout.title = "FID"

    # ms
    ms.processor.add(
        ca.p.edit.ReplaceSpans(value=0,
                               x_spans=(
                                   (8.9, 10),  # decane
                                   (44.1, 45.1)  # PPh3
                               ), invert=True),
        ca.p.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4)
    )
    peak_locations = ca.a.peak_picking.find_peaks_scipy(ms, scipy_kwargs={"height": 8000, "width": 0.1})
    peaks = ca.a.integration.rolling_ball(peak_locations, n=5, min_height=0.002, n_points_with_pos_slope=2)
    ms_compounds = ca.a.ms_analysis.search_by_retention_time(picking_lib_ms, peaks)
    
    ms_fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    ca.plotting.signal(ms, fig=ms_fig)
    ca.plotting.peaks(ms_compounds, fig=ms_fig)
    ms_fig.layout.title = "MS"

    print("finished analyzing:", label)
    return ms_fig, fid_fig, ms_compounds, fid_compounds


def get_figure_path(root_folder: pathlib.Path) -> pathlib.Path:
    import shutil
    figure_folder = root_folder / "figs"
    if figure_folder.exists():
        shutil.rmtree(figure_folder)
    figure_folder.mkdir(parents=True, exist_ok=True)

    return figure_folder


def process_timeseries(data_path: pathlib.Path, labels: list[str], times: np.ndarray):
    figure_folder = get_figure_path(data_path)

    # process data
    fid_compounds, ms_compounds, fid_figs, ms_figs = [], [], [], []
    for label in labels:
        ms_fig_, fid_fig_, ms_comp, fid_comp = process_single(data_path, label)
        ms_figs.append(ms_fig_)
        fid_figs.append(fid_fig_)
        ms_compounds.append(ms_comp)
        fid_compounds.append(fid_comp)

    # process timeseries
    fid_timeseries = ResultTimeSeries("decane_fid", INTERNAL_STANDARD, internal_standard_mmol=0.0109)
    ms_timeseries = ResultTimeSeries("decane_ms", INTERNAL_STANDARD, internal_standard_mmol=0.0109)
    for i in range(len(times)):
        fid_timeseries.add_result(fid_compounds[i], times[i])
        ms_timeseries.add_result(ms_compounds[i], times[i])


    decane_conc = 0.092

    fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    plot_results(fid_timeseries, PLOTTING_GROUPS, fig=fig, add_zero=True)
    fig.add_scatter(x=[0, times[-1]], y=[decane_conc, decane_conc], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    fig.layout.xaxis.title = "<b>time (min)<br>"
    fig.layout.yaxis.title = "<b>mmol<br>"
    fid_figs.append(fig)
    fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    plot_results(ms_timeseries, PLOTTING_GROUPS, fig=fig, add_zero=True)
    fig.add_scatter(x=[0, times[-1]], y=[decane_conc, decane_conc], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    fig.layout.xaxis.title = "<b>time (min)<br>"
    fig.layout.yaxis.title = "<b>mmol<br>"
    ms_figs.append(fig)

    fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    plot_results(fid_timeseries, PLOTTING_GROUPS, fig=fig, carbon=True)
    fig.add_scatter(x=[0, times[-1]], y=[decane_conc*10, decane_conc*10], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    fig.layout.xaxis.title = "<b>time (min)<br>"
    fig.layout.yaxis.title = "<b>mmol of carbon<br>"
    fid_figs.append(fig)
    fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    plot_results(ms_timeseries, PLOTTING_GROUPS, fig=fig, carbon=True)
    fig.add_scatter(x=[0, times[-1]], y=[decane_conc*10, decane_conc*10], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    fig.layout.xaxis.title = "<b>time (min)<br>"
    fig.layout.yaxis.title = "<b>mmol of carbon<br>"
    ms_figs.append(fig)

    ca.plotting.plotly_utils.merge_figures(fid_figs, filename=figure_folder / "fid")
    ca.plotting.plotly_utils.merge_figures(ms_figs, filename=figure_folder / "ms")

    data = ms_timeseries.to_csv_str()
    print(data)


def main():
    data_path = pathlib.Path(r"C:\Users\nicep\Desktop\10_10_reduce")
    times = np.arange(1, 12)
    labels = [f"DJW-10-10-{i}h-TMS" for i in times]
    process_timeseries(data_path, labels, times)


if __name__ == "__main__":
    main()
