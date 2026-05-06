import pathlib

import numpy as np
import polars as pl
import plotly.graph_objs as go

import chem_analysis as ca

import development.plotly_utils as plotly_utils

folder_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\publications\DW_AA_oxidation\data\hydro_peroxide_N2\GC-raw")
data_label = "DJW-11-98-t"
times = np.array([0, 15, 30, 60, 120, 240])


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    return ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)


labels_TMS = [data_label + str(i) + "-TMS" for i in times]
labels_PPh = [data_label + str(i) + "-PPh3" for i in times]
raw_data_TMS: list[tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]] = [load_single_file(folder_path, label) for label in labels_TMS]
raw_data_PPh: list[tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]] = [load_single_file(folder_path, label) for label in labels_PPh]


def plot_ms_fid(i: int, ms: ca.gc_lc.GCMSSignal, fid: ca.gc_lc.GCSignal, peaks: list[ca.a.PeakData]):
    fig_ms = ca.plotting.signal(ms)

    fig_fid = ca.plotting.signal(fid)
    fig_fid = ca.plotting.add_peaks(peaks, fig=fig_fid, colors="multi")
    fig_fid.layout.yaxis.range = [-50_000, 4_000_000]

    plotly_utils.merge_html_figs([fig_ms, fig_fid], filename=f"plots_{i}_{times[i]}.html", auto_open=True)


def tn_TMS(i: int):
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data_TMS[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [29.3, 30], label="TCB", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [25.1, 25.7], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.91, 26.4], label="A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.5], label="AA10_3", to_zero=True)
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def tn_PPh(i: int):
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data_PPh[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [29.3, 30], label="TCB_PPh", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.1, 25.7], label="K10_3_PPh", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.91, 26.4], label="A10_3_PPh", to_zero=True),

    ]

    # plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def plot_series():
    fig = go.Figure()
    colors = plotly_utils.get_colors_from_scale(len(raw_data_TMS), "Viridis")
    colors.reverse()
    for i, data in enumerate(raw_data_PPh):
        x, y = data[1].x, data[1].y_normalized_by_max([28, 30])
        y = y  # + 0.02*i
        fig.add_scatter(x=x, y=y, mode="lines", line_color=colors[i], name=f"{times[i]} min")

    # full
    # fig.layout.xaxis.range = 5.5, None
    # fig.layout.yaxis.range = -0.05, 2

    # zoom 1
    fig.layout.xaxis.range = 5, 35
    fig.layout.yaxis.range = -0.02, 1.2

    fig.layout.xaxis.title.text = "retention time (min)"
    fig.layout.yaxis.title.text = "FID intensity"
    fig.layout.margin.t = 20
    fig.show('browser')
    # fig.write_image("DJW-11-98-PPh.png")  # _zoom1


def main():
    data = [tn_TMS(i) for i in range(len(times))]
    data2 = [tn_PPh(i) for i in range(len(times))]
    for i in range(len(times)):
        data[i] |= data2[i]

    df = pl.DataFrame(data)
    df.with_columns(pl.col(col).fill_null(0) for col in df.columns)
    print(df)
    df.write_csv(data_label + "_areas.csv")

    # print(t6())


if __name__ == "__main__":
    # plot_series()
    main()
