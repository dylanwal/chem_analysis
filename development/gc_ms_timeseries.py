import glob
import pathlib
import re
import shutil

import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json"
LIBRARY = ca.mass_spec.GCLibrary.from_JSON(lib_path)
INTERNAL_STANDARD = LIBRARY.find_by_label("TCB")
picking_lib_fid = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_fid")
picking_lib_ms = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_ms")


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


def peaks_to_compounds(
        results: list[ca.analysis.peak_picking.ResultPeaks],
        picking_library: ca.analysis.ms_analysis.PickingLibrary
):
    compounds = []
    for result in results:
        compounds.append(ca.analysis.ms_analysis.search_retention_time(picking_library, result))

    return compounds


def get_figure_path(root_folder: pathlib.Path) -> pathlib.Path:
    figure_folder = root_folder / "figs"
    if figure_folder.exists():
        shutil.rmtree(figure_folder)
    figure_folder.mkdir(parents=True, exist_ok=True)

    return figure_folder


def get_data(data_path: pathlib.Path, pattern: str, sort_func=None):
    if sort_func is None:
        pattern_compile = re.compile(pattern.replace("*", "([0-9]+)"))
        sort_func = lambda x: int(pattern_compile.match(x).groups()[0])

    specific_folders = glob.glob(pattern, root_dir=data_path)
    print(len(specific_folders), "files found for analysis")
    specific_folders.sort(key=sort_func)
    specific_folders = [data_path / file for file in specific_folders]
    return [ca.gc_lc.GCParser.from_Agilent_D_folder(folder) for folder in specific_folders]


def process_timeseries(data_path: str, pattern: str, peak_mask):
    if isinstance(data_path, str):
        data_path = pathlib.Path(data_path)
    figure_folder = get_figure_path(data_path)
    data = get_data(data_path, pattern)

    # process data
    fid_peak_results = []
    ms_peak_results = []
    for ms, fid in data:
        fid_peak_results.append(process_one_signal(fid, "fid", figure_folder, peak_mask))
        ms_peak_results.append(process_one_signal(ms, "ms", figure_folder, peak_mask))

    # peaks --> compounds
    fid_compounds = peaks_to_compounds(fid_peak_results, picking_lib_fid)
    ms_compounds = peaks_to_compounds(ms_peak_results, picking_lib_ms)

    # plotting timeseries



def main():
    data_path = r"C:\Users\nicep\Desktop\10_13"
    pattern = "DJW-10-13-*min-PPh3.D"
    peak_mask = ca.processing.weigths.Spans([4, 31], invert=True)

    process_timeseries(data_path, pattern, peak_mask)


if __name__ == "__main__":
    main()
