import pathlib

import numpy as np

import chem_analysis as ca

from result_series import results_to_timeseries

lib_FID = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\intergate_FID_AA.txt"
parameters = ""


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    return ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)


def process_fid(sig: ca.gc_lc.GCSignal):
    proc = ca.p.Processor(
        ca.p.baseline.MorphologicalAverage()
    )
    proc_sig = proc.run(sig)

    peaks = chem_analysis.analysis.peaks.integration.integrate_from_file(proc_sig, lib_FID)

    fig = ca.plot.signal(proc_sig)
    fig = ca.plot.peaks(peaks, fig=fig)
    max_y = max(peak.properties.max_y for peak in peaks.peaks)
    fig.layout.yaxis.range = [-0.1 * max_y, 1.1 * max_y]
    fig.layout.title = "FID"
    return fig, peaks


def process_ms(sig: ca.gc_lc.GCMSSignal):
    fig = ca.plot.signal(sig)
    fig.layout.title = "MS"
    return fig


def process_timeseries(data_path: pathlib.Path, data_label: str, labels: list[str], times: np.ndarray):
    global parameters
    parameters += str(data_path) + "\n" + str(data_label) + "\n"

    # process data
    fid_compounds, fid_figs, ms_figs = [], [], []
    for label in labels:
        ms, fid = load_single_file(data_path / "GCMS", label)
        fid_fig_, fid_comp = process_fid(fid)
        fid_figs.append(fid_fig_)
        fid_compounds.append(fid_comp)
        ms_figs.append(process_ms(ms))
        print(f"processed: {label}")

    # re-organizing data
    fid_timeseries = results_to_timeseries(fid_compounds, times)

    # saving data
    ca.plotting.plotly_utils.merge_figures(fid_figs, filename=data_path / (data_label + '_fid.html'))
    ca.plotting.plotly_utils.merge_figures(ms_figs, filename=data_path / (data_label + '_ms.html'))
    data = fid_timeseries.write_csv(data_path / (data_label + '_fid.csv'))
    print(data)


def main_baseline():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_65")
    times = np.array([30, 60, 120, 240, 360, 720, 1950])  # 30, 60, 120, 240, 360, 720, 1950
    data_label = "DJW-11-65-v1"
    labels = [data_label + f"-t{i}-TMS" for i in times]

    figs = []
    for label in labels:
        ms, fid = load_single_file(data_path / "GCMS", label)
        proc = ca.p.Processor(
            ca.p.baseline.MorphologicalAverage()
        )
        proc_sig = proc.run(fid)
        # print(f"window: {proc.methods[0].window_size}")
        peaks = chem_analysis.analysis.peaks.integration.integrate_from_file(proc_sig, lib_FID)
        fig = ca.plot.signal(fid)
        fig = ca.plot.baseline(proc, fig=fig)
        # fig = ca.plot.peaks(peaks, fig=fig)
        max_y = max(peak.properties.max_y for peak in peaks.peaks)
        fig.layout.yaxis.range = [-0.1 * max_y, 1.1 * max_y]
        fig.layout.title = "FID"

        figs.append(fig)
        print(f"processed: {label}")

    ca.plotting.plotly_utils.merge_figures(figs, filename=data_path / 'baseline_fid.html')


def main():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_65")
    times = np.array([30, 60, 120, 240, 360, 720, 1950])  # 30, 60, 120, 240, 360, 720, 1950
    data_label = "DJW-11-65-v1"
    labels = [data_label + f"-t{i}-TMS" for i in times]
    process_timeseries(data_path, data_label, labels, times)


def main_plot():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_65")
    times = np.array([30, 60, 120, 240, 360, 720, 1950])  # 30, 60, 120, 240, 360, 720, 1950
    data_label = "DJW-11-65-v1"
    labels = [data_label + f"-t{i}-TMS" for i in times]
    sigs = [load_single_file(data_path / "GCMS", label) for label in labels]
    fids = ca.gc_lc.GCSignal2D.from_signals([i[1] for i in sigs], times, "time",
                                            unify_method=ca.base_obj.UnifyMethodExpandInterpolate())
    mss = [i[0] for i in sigs]

    fig = ca.plot.signal2D_overlap_signals(fids)
    fig.show()


if __name__ == "__main__":
    main_baseline()
    # main()
    # main_plot()
