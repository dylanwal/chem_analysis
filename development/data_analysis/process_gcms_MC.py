import pathlib

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca

from result_series import results_to_timeseries

lib_FID = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\intergate_FID_MC.txt"
parameters = ""


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

    peaks_fid = ca.a.integration.integrate_from_file(fid, lib_FID)
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
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_53")
    times = np.array([0, 15, 30, 60, 120, 240])
    data_label = "DJW-11-53-run1"
    labels = [data_label + f"-t{i}" for i in times]
    process_timeseries(data_path, data_label, labels, times)


if __name__ == "__main__":
    main()
