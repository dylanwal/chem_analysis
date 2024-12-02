import pathlib

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca

from result_series import results_to_timeseries

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library_color.json"
LIBRARY = ca.library.Library.from_JSON(lib_path)
picking_lib_fid = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_fid")
picking_lib_ms = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_ms")

parameters = ""
lib_FID = r"C:\Users\nicep\Desktop\intergate_FID_MC.txt"


def process_single(data_path: pathlib.Path, label: str):
    # load data
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)

    # fid
    baseline_proc = ca.p.baseline.BaselineWithMask(
        baseline_method=ca.p.baseline.Polynomial(degree=0),
        mask=ca.p.weigths.Spans((None, 1.5))
    )
    fid.processor.add(baseline_proc)

    # peak_locations = ca.a.peak_picking.find_peaks_scipy(fid,
    #                                                     mask=ca.p.weigths.Spans((4.1, 44)),
    #                                                     scipy_kwargs={"height": 20_000, "width": 0.15}
    #                                                     )
    # peaks = ca.a.integration.integrate_by_fitting_normal_distribution(peak_locations)
    # peaks = ca.a.integration.rolling_ball(peak_locations, n=5, min_height=0.002, n_points_with_pos_slope=2)
    peaks_fid = ca.a.integration.integrate_from_file_trapz(fid, lib_FID)
    # fid_compounds = ca.a.ms_analysis.search_by_retention_time(picking_lib_fid, peaks)

    fid_fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    ca.plotting.signal(fid, fig=fid_fig)
    ca.plotting.peaks(peaks_fid, fig=fid_fig)
    max_y = max(peak.properties.max_y for peak in peaks_fid.peaks)
    fid_fig.layout.yaxis.range = [-0.1*max_y, 1.1*max_y]
    fid_fig.layout.title = "FID"

    print("finished analyzing:", label)
    global parameters
    parameters += "\n" + str(data_path) + "\n\t fid:" + str(fid.processor.methods)
    return fid_fig, peaks_fid


def process_timeseries(data_path: pathlib.Path, data_label: str, labels: list[str], times: np.ndarray):
    global parameters
    parameters += str(data_path) + "\n" + str(data_label) + "\n"

    # process data
    fid_compounds, ms_compounds, fid_figs, ms_figs = [], [], [], []
    for label in labels:
        fid_fig_, fid_comp = process_single(data_path / "GCMS", label)
        fid_figs.append(fid_fig_)
        fid_compounds.append(fid_comp)

    # re-organizing data
    fid_timeseries = results_to_timeseries(fid_compounds, times)

    # saving data
    ca.plotting.plotly_utils.merge_figures(fid_figs, filename=data_path / (data_label + '_fid.html'))
    ca.plotting.plotly_utils.merge_figures(ms_figs, filename=data_path / (data_label + '_ms.html'))
    data = fid_timeseries.to_csv(data_path / (data_label + '_fid.csv'))

    with open(data_path.parent / (data_label + "_params.txt"), mode='w') as f:
        f.write(parameters)

    print(data)


def main():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_65")
    times = np.array([30, 60, 120, 240, 360, 720, 1950])  # 30, 60, 120, 240, 360, 720, 1950
    data_label = "DJW-11-65-v1"
    labels = [data_label + f"-t{i}-TMS" for i in times]
    process_timeseries(data_path, data_label, labels, times)


if __name__ == "__main__":
    main()
