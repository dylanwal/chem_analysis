import glob
import pathlib
import shutil

import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json"
LIBRARY = ca.mass_spec.GCLibrary.from_JSON(lib_path)


def process_one_signal(signal: ca.base_obj.Signal, type_: str, figure_folder, peak_mask) \
        -> ca.analysis.peak_picking.ResultPeaks:
    if type_ == "fid":
        signal.processor.add(ca.processing.edit.ReplaceSpans(value=0, x_spans=(1.3, 3.7), invert=True))
    signal.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True))

    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(signal,
                                                               mask=peak_mask,
                                                               scipy_kwargs={"height": 6000, "width": 0.1}
                                                               )
    peak_result = ca.analysis.boundary_detection.rolling_ball(peak_locations, n=5, min_height=0.002,
                                                              n_points_with_pos_slope=2)

    # plotting peak results
    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    # ca.plotting.baseline(signal, fig=fig)
    ca.plotting.signal(signal, fig=fig)
    ca.plotting.peaks(peak_result, fig=fig)

    fig.write_html(figure_folder / f"signal_{type_}_{signal.name}.html", include_plotlyjs='cdn')

    return peak_result


def get_figure_path(root_folder: str) -> pathlib.Path:
    figure_folder = pathlib.Path(root_folder) / "figs"
    if figure_folder.exists():
        shutil.rmtree(figure_folder)
    figure_folder.mkdir(parents=True, exist_ok=True)

    return figure_folder


def get_data(data_path, sort_func):
    specific_folders = glob.glob(data_path + "/*.D")
    print(len(specific_folders), "files found for analysis")
    specific_folders.sort(key=sort_func)
    return [ca.gc_lc.GCParser.from_Agilent_D_folder(folder) for folder in specific_folders]


def process_timeseries(data_path: str, sort_func, peak_mask):
    figure_folder = get_figure_path(data_path)
    data = get_data(data_path, sort_func)

    # process data
    fid_peak_results = []
    ms_peak_results = []
    for ms, fid in data:
        fid_peak_results.append(process_one_signal(fid, "fid", figure_folder, peak_mask))
        ms_peak_results.append(process_one_signal(ms, "ms", figure_folder, peak_mask))

    # peaks --> compounds
    picking_lib = LIBRARY.to_picking_library()
    peaks_compounds = ca.analysis.peak_picking.library_search.find_peaks_retention_time_library(ms_peak_results,
                                                                                                picking_lib)
    # plotting timeseries


def main():
    data_path = ""
    sort_func = lambda x: x.split("//")[-1].replace("redo", "").replace("-", "").replace("_", "")
    peak_mask = ca.processing.weigths.Spans([4, 31], invert=True)
    process_timeseries(data_path, sort_func, peak_mask)


if __name__ == "__main__":
    main()
