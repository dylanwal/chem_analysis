import pathlib

import numpy as np
import polars as pl
import plotly.graph_objs as go

import chem_analysis as ca

import development.plotly_utils as plotly_utils

folder_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\publications\DW_AA_oxidation\data\ketone_stablity\GC-raw")
data_label = "ZS-7-V2-t"
times = np.array([0, 30, 60, 120, 240, 480, 1440])  #


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    return ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)


labels = [data_label + str(i) for i in times]
raw_data: list[tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]] = [load_single_file(folder_path, label) for label in labels]


def plot_ms_fid(i: int, ms: ca.gc_lc.GCMSSignal, fid: ca.gc_lc.GCSignal, peaks: list[ca.a.PeakData]):
    fig_ms = ca.plotting.signal(ms)

    fig_fid = ca.plotting.signal(fid)
    fig_fid = ca.plotting.add_peaks(peaks, fig=fig_fid, colors="multi")
    fig_fid.layout.yaxis.range = [-50_000, 4_000_000]

    plotly_utils.merge_html_figs([fig_ms, fig_fid], filename=f"plots_{i}_{times[i]}.html", auto_open=True)


def tn(i: int):
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [29.03, 29.71], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.58, 8.9], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.45, 25.8], label="K10_2", to_zero=True)
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def plot_series():
    fig = go.Figure()
    colors = plotly_utils.get_colors_from_scale(len(raw_data), "Viridis")
    colors.reverse()
    for i, data in enumerate(raw_data):
        x, y = data[1].x, data[1].y_normalized_by_max([28, 30])
        y = y  # + 0.02*i
        fig.add_scatter(x=x, y=y, mode="lines", line_color=colors[i], name=f"{times[i]} min")

    # full
    fig.layout.xaxis.range = 5.5, None
    fig.layout.yaxis.range = -0.05, 2

    # zoom 1
    fig.layout.xaxis.range = 5.5, 30
    fig.layout.yaxis.range = -0.02, 2

    fig.layout.xaxis.title.text = "retention time (min)"
    fig.layout.yaxis.title.text = "FID intensity"
    fig.layout.margin.t = 20
    fig.show('browser')
    # fig.write_image("ZS-7-V2.png")  # _zoom1


def main():
    data = [tn(i) for i in range(len(times))]
    df = pl.DataFrame(data)
    df.with_columns(pl.col(col).fill_null(0) for col in df.columns)
    print(df)
    df.write_csv(data_label + "_areas.csv")

    # print(t6())


if __name__ == "__main__":
    # plot_series()
    main()
