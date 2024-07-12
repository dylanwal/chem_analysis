import pathlib
import re
import shutil

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


def process_single(data_path):
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_folder(data_path)

    # fid
    fid.processor.add(ca.processing.edit.ReplaceSpans(value=0, x_spans=(1.3, 3.7), invert=True))
    fid.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
    )
    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(fid,
                                                               scipy_kwargs={"height": 8000, "width": 0.1}
                                                               )
    fid_peaks = ca.analysis.integration.rolling_ball(peak_locations, n=5, min_height=0.002,
                                                     n_points_with_pos_slope=2)
    fid_compounds = ca.analysis.ms_analysis.search_by_retention_time(picking_lib_fid, fid_peaks)

    fid_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(fid, fig=fid_fig)
    ca.plotting.peaks(fid_compounds, fig=fid_fig)
    fid_fig.layout.title = "FID"

    # ms
    ms.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
    )
    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(ms,
                                                               scipy_kwargs={"prominence": 10000}
                                                               )
    ms_peaks = ca.analysis.integration.rolling_ball(peak_locations, n=5, min_height=0.002,
                                                    n_points_with_pos_slope=2)
    ms_compounds = ca.analysis.ms_analysis.search_by_retention_time(picking_lib_ms, ms_peaks)

    # plotting peak results
    ms_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(ms, fig=ms_fig)
    ca.plotting.peaks(ms_compounds, fig=ms_fig)
    ms_fig.layout.title = "MS"

    print("finished analyzing:", data_path.name)
    return ms_fig, fid_fig, ms_compounds, fid_compounds


def get_figure_path(root_folder: pathlib.Path) -> pathlib.Path:
    figure_folder = root_folder / "figs"
    if figure_folder.exists():
        shutil.rmtree(figure_folder)
    figure_folder.mkdir(parents=True, exist_ok=True)

    return figure_folder


def get_folders(data_path: pathlib.Path, pattern: str) -> tuple[list[pathlib.Path], np.ndarray]:
    pattern_compile = re.compile(pattern.replace("*", "([0-9]+)"))

    def sort_func(x: pathlib.Path) -> int:
        return int(pattern_compile.match(x.name).groups()[0])

    specific_folders = list(data_path.glob(pattern))
    print(len(specific_folders), "files found for analysis")

    specific_folders.sort(key=sort_func)
    times_ = [sort_func(folder) for folder in specific_folders]

    return [data_path / file for file in specific_folders], np.array(times_)


def process_timeseries(data_path: str, pattern: str):
    if isinstance(data_path, str):
        data_path = pathlib.Path(data_path)
    figure_folder = get_figure_path(data_path)
    folders, times = get_folders(data_path, pattern)

    # process data
    fid_compounds, ms_compounds, fid_figs, ms_figs = [], [], [], []
    for folder in folders:
        ms_fig_, fid_fig_, ms_comp, fid_comp = process_single(folder)
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

    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    plot_results(fid_timeseries, PLOTTING_GROUPS, fig=fig, add_zero=True)
    fig.add_scatter(x=[0, 20], y=[0.0658, 0.0658], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    fig.layout.xaxis.title = "<b>time (min)<br>"
    fig.layout.yaxis.title = "<b>mmol<br>"
    fid_figs.append(fig)
    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    plot_results(ms_timeseries, PLOTTING_GROUPS, fig=fig, add_zero=True)
    fig.add_scatter(x=[0, 20], y=[0.0658, 0.0658], mode="lines", line={"color": "black", "dash": "dash"}, name="decane_init")
    fig.layout.xaxis.title = "<b>time (min)<br>"
    fig.layout.yaxis.title = "<b>mmol<br>"
    ms_figs.append(fig)

    ca.plotting.PlotlyConfig.merge_figures(fid_figs, filename=figure_folder / "fid")
    ca.plotting.PlotlyConfig.merge_figures(ms_figs, filename=figure_folder / "ms")

    data = ms_timeseries.to_csv_str()
    print(data)


def main():
    data_path = r"C:\Users\nicep\Desktop\11_22"
    pattern = "DJW-11-22-*min-TMS.D"
    process_timeseries(data_path, pattern)


if __name__ == "__main__":
    main()
