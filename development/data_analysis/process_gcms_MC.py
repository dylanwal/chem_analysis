import pathlib

import numpy as np
import plotly.graph_objs as go
from scipy.ndimage import gaussian_filter1d

import chem_analysis as ca

from result_series import results_to_timeseries

lib_FID = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\intergate_FID_MC.txt"
parameters = ""


def build_baseline(data_path: pathlib.Path, label: str):
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)

    baseline = ca.p.baseline.AdaptiveAsymmetricLeastSquared(save_result=True)
    fid.processor.add(
        ca.p.edit.CutSpans(((1.2, 11), (29, 31)), invert=True),
        ca.p.resampling.EveryN(step=[int(fid.x_raw.size/6000)]),
        ca.p.smoothing.Gaussian(sigma=10),
        baseline
    )

    fig = ca.plot.signal(fid, raw=True)
    fig = ca.plot.baseline(fid, fig=fig)
    min_, max_ = np.min(baseline.baseline), np.max(baseline.baseline)
    fig.layout.yaxis.range = min_ - (max_ - min_) / 4, max_ + (max_ - min_) / 2

    return baseline.x, baseline.baseline, fig


def process_single(data_path: pathlib.Path, label: str, x_baseline: np.ndarray, y_baseline: np.ndarray):
    # load data
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)

    # fid
    # baseline_proc = ca.p.baseline.BaselineWithMask(
    #     baseline_method=ca.p.baseline.Polynomial(degree=0),
    #     mask=ca.p.weigths.Spans((None, 1.5))
    # )
    # fid.processor.add(baseline_proc)
    fid.processor.add(ca.p.translations.AlignMaxValue((29.5, 29.7), x_value=29.566))
    fid.processor.add(ca.p.baseline.Subtract(y=y_baseline, x=x_baseline))
    fid.processor.add(ca.p.baseline.BaselineWithMask(
        ca.p.baseline.AdaptiveAsymmetricLeastSquared(lambda_=1e2),
        mask=ca.p.weigths.SlidingWindowStd(window=10, sections=200, smoother=lambda x: gaussian_filter1d(x, 5))
    ))

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


def process_timeseries(data_path: pathlib.Path, data_label: str, labels: list[str], times: np.ndarray, multiple: bool = False):
    global parameters
    parameters += str(data_path) + "\n" + str(data_label) + "\n"
    figs = []

    # get baseline
    x, y, fig = build_baseline(data_path / "GCMS", labels[0])
    figs.append(fig)

    # process data
    fid_compounds = []
    for label in labels:
        fid_fig, fid_comp = process_single(data_path / "GCMS", label, x_baseline=x, y_baseline=y)
        figs.append(fid_fig)
        fid_compounds.append(fid_comp)

    # re-organizing data
    fid_timeseries = results_to_timeseries(fid_compounds, times)

    # saving data
    if multiple:
        ca.plotting.plotly_utils.merge_figures(figs, filename=data_path / "figs" / (data_label + '_figs.html'))
        data = fid_timeseries.to_csv(data_path / "areas" / (data_label + '_fid.csv'))
    else:
        ca.plotting.plotly_utils.merge_figures(figs, filename=data_path / (data_label + '_figs.html'))
        data = fid_timeseries.to_csv(data_path / (data_label + '_fid.csv'))

    with open(data_path.parent / (data_label + "_params.txt"), mode='w') as f:
        f.write(parameters)

    print(data)


def main_multi():
    data_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\11\11_53")

    # get all data_labels and times
    import glob
    import re
    files = glob.glob(str(data_path / "GCMS\*.ms"))
    labels = dict()
    for file_name in files:
        integers = re.findall(r'\d+', file_name)
        integers = integers[slice(5, None)]
        run_num, time_ = [int(i) for i in integers]
        if run_num in labels:
            labels[run_num].append(time_)
        else:
            labels[run_num] = [time_]

    # run_xy method on everything
    import multiprocessing
    tasks = [(run_num, times, data_path) for run_num, times in labels.items()]
    with multiprocessing.Pool(multiprocessing.cpu_count()-1) as pool:
        pool.map(process_run, tasks)


def process_run(args):
    run_num, times, data_path = args
    times = np.array(times)
    times.sort()
    data_label = f"DJW-11-53-run_xy{run_num}"
    labels = [data_label + f"-t{i}" for i in times]
    process_timeseries(data_path, data_label, labels, times, True)
    print(f"Finished analyzing run_xy {run_num}")


def main():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_53")
    times = np.array([0, 15, 30, 60, 120, 240])
    data_label = "DJW-11-53-run1"
    labels = [data_label + f"-t{i}" for i in times]
    process_timeseries(data_path, data_label, labels, times)


if __name__ == "__main__":
    # main()
    main_multi()
