import pathlib

import numpy as np
import polars as pl
import plotly.graph_objs as go

import chem_analysis as ca

import development.plotly_utils as plotly_utils

folder_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\11\11_115\GCMS")
data_label = "DJW-11-115-V1-t"
times = np.array([0, 30, 60, 90, 120, 150, 180, 240, 300, 360, 480, 600])  #


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    return ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)


labels = [data_label + str(i) for i in times]
raw_data: list[tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]] = [load_single_file(folder_path, label) for label in
                                                                 labels if "600" not in label]


def load_600() -> tuple[ca.gc_lc.GCSignal, ca.gc_lc.GCSignal]:
    ms_data = np.loadtxt(r"C:\Users\nicep\Desktop\research_wis\data\11\11_115\GCMS\DJW-11-115-V1-t600_ms.csv",
                         delimiter=",", skiprows=2)
    ms = ca.gc_lc.GCSignal(x=ms_data[:, 1], y=ms_data[:, 2])

    fid_data = np.loadtxt(r"C:\Users\nicep\Desktop\research_wis\data\11\11_115\GCMS\DJW-11-115-V1-t600_fid.csv",
                          delimiter=",", skiprows=2)
    fid = ca.gc_lc.GCSignal(x=fid_data[:, 1], y=fid_data[:, 2])

    return ms, fid


raw_data.insert(11, load_600())


def plot_ms_fid(i: int, ms: ca.gc_lc.GCMSSignal, fid: ca.gc_lc.GCSignal, peaks: list[ca.a.PeakData]):
    fig_ms = ca.plotting.signal(ms)

    fig_fid = ca.plotting.signal(fid)
    fig_fid = ca.plotting.add_peaks(peaks, fig=fig_fid, colors="multi")
    fig_fid.layout.yaxis.range = [-50_000, 4_000_000]

    plotly_utils.merge_html_figs([fig_ms, fig_fid], filename=f"plots_{i}_{times[i]}.html", auto_open=True)


def t0():
    i = 0
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [28.96, 29.64], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True)
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t1():
    i = 1
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [28.96, 29.27], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True)
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t2():
    i = 2
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [6.2, 6.40], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True),
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t3():
    i = 3
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.97, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.4, 30.6], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True), ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t4():
    i = 4
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.97, 29.29], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.29, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.37, 30.57], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True), ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t5():
    i = 5
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.95, 29.28], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.28, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True), ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t6():
    i = 6
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.95, 29.285], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.285, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True), ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t7():
    i = 7
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.94, 29.285], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.285, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True), ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t8():
    i = 8
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.94, 29.28], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.28, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True), ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t9():
    i = 9
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.94, 29.283], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.283, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True), ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t10():
    i = 10
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.94, 29.28], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.28, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True), ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t11():
    i = 11
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks = [
        ca.a.peaks.peak_from_span(fid, [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.30, 8.48], label="K7_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [28.94, 29.28], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.28, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.46, 43.68], label="DA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [46.5, 46.72], label="DA7", to_zero=True),
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


# def t12():
#     i=12
#     data = dict()
#     data["time"] = times[i]
#     ms, fid = raw_data[i]
#     peaks =  [
#         ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
#         ca.a.peaks.peak_from_span(fid,  [28.94, 29.28], label="TCB", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [29.28, 29.51], label="KA6_3", to_zero=True),
#         # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA6", to_zero=True),  #mix
#     ]
#
#     plot_ms_fid(i, ms, fid, peaks)
#     data.update({peak.label: peak.properties.area for peak in peaks})
#     return data
#
# def t13():
#     i=13
#     data = dict()
#     data["time"] = times[i]
#     ms, fid = raw_data[i]
#     peaks =  [
#         ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
#         ca.a.peaks.peak_from_span(fid,  [28.94, 29.27], label="TCB", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [29.27, 29.51], label="KA6_3", to_zero=True),
#         # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA6", to_zero=True),  #mix
#     ]
#
#     plot_ms_fid(i, ms, fid, peaks)
#     data.update({peak.label: peak.properties.area for peak in peaks})
#     return data
#
# def t14():
#     i=14
#     data = dict()
#     data["time"] = times[i]
#     ms, fid = raw_data[i]
#     peaks =  [
#         ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [16.6, 16.9], label="L5_4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [17.6, 17.9], label="HAA2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.22, 19.4], label="K9_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.44, 19.63], label="HA3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [19.63, 20], label="HA3_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A10_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
#         ca.a.peaks.peak_from_span(fid,  [28.94, 29.25], label="TCB", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [29.25, 29.51], label="KA6_3", to_zero=True),
#         # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
#         ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA6", to_zero=True),  #mix
#     ]
#
#     plot_ms_fid(i, ms, fid, peaks)
#     data.update({peak.label: peak.properties.area for peak in peaks})
#     return data


KA5_2_ms_raw = [
    ['36.9', '4.84'],
['37.1', '193.29'],
['38.1', '530.04'],
['39.1', '3728.16'],
['40.1', '724.27'],
['41.1', '3615.08'],
['42.1', '11009.32'],
['43.1', '130223.66'],
['44.1', '11714.63'],
['45.1', '67817.89'],
['46.1', '5130.41'],
['47.1', '35377.57'],
['48.', '1837.18'],
['49.', '1379.36'],
['50.1', '648.31'],
['51.1', '5767.51'],
['52.', '786.26'],
['53.1', '11994.91'],
['54.1', '2332.98'],
['55.1', '34358.42'],
['56.1', '29994.85'],
['57.1', '8642.11'],
['58.1', '18245.63'],
['59.1', '19887.65'],
['60.1', '8562.45'],
['61.1', '37027.8'],
['62.', '2767.32'],
['63.', '1736.46'],
['64.', '67.48'],
['64.2', '5.91'],
['65.1', '1410.23'],
['66.', '358.87'],
['67.', '569.22'],
['68.1', '336.64'],
['69.1', '2002.15'],
['70.1', '4257.57'],
['71.1', '13907.4'],
['72.1', '12121.98'],
['73.1', '304192.06'],
['74.1', '39577.13'],
['75.1', '672143.13'],
['76.1', '52099.11'],
['77.1', '35565.32'],
['78.1', '1360.29'],
['79.1', '4104.53'],
['79.9', '158.22'],
['80.', '23.64'],
['81.1', '24440.27'],
['82.1', '1671.87'],
['83.', '1565.67'],
['84.1', '574.36'],
['85.1', '5966.47'],
['86.1', '20228.15'],
['87.1', '5600.48'],
['88.1', '2144.44'],
['89.1', '3168.35'],
['90.1', '433.27'],
['91.1', '6413.84'],
['92.1', '581.'],
['93.1', '1554.8'],
['94.', '109.14'],
['95.', '751.16'],
['96.', '221.47'],
['97.1', '1356.54'],
['98.1', '19777.36'],
['99.1', '45124.21'],
['100.1', '3237.51'],
['101.1', '5511.53'],
['102.1', '1474.69'],
['103.1', '3297.6'],
['104.1', '342.26'],
['105.1', '695.13'],
['106.1', '5.04'],
['109.1', '1146.07'],
['110.2', '12.06'],
['111.1', '22899.33'],
['112.1', '3150.96'],
['113.1', '12415.29'],
['114.1', '1577.12'],
['115.1', '7047.43'],
['116.1', '3017.59'],
['117.1', '7963.84'],
['118.', '1104.48'],
['119.', '473.46'],
['125.1', '541.04'],
['126.1', '20.23'],
['127.1', '12518.62'],
['128.1', '2644.15'],
['129.1', '66628.32'],
['130.1', '32310.53'],
['131.1', '72915.98'],
['132.1', '9044.39'],
['133.1', '3910.65'],
['134.1', '291.59'],
['134.1', '38.48'],
['135.', '39.04'],
['137.1', '9.93'],
['139.', '5.41'],
['141.1', '51.39'],
['142.2', '12.78'],
['143.1', '2301.1'],
['144.2', '18425.5'],
['145.1', '367965.'],
['146.1', '108845.14'],
['147.1', '25124.06'],
['148.1', '4211.28'],
['149.1', '332.24'],
['155.1', '39869.84'],
['156.1', '5356.64'],
['157.1', '2120.24'],
['158.1', '248.37'],
['159.', '4.86'],
['160.1', '303.2'],
['170.1', '2909.98'],
['171.1', '609.43'],
['173.1', '400330.16'],
['174.1', '53958.93'],
['175.1', '19674.31'],
['176.1', '1600.07'],
['177.1', '105.53'],
['177.1', '59.73'],
['188.1', '7570.45'],
['189.1', '1130.02'],
['190.1', '398.54'],
['191.', '84.45'],
]

# getting the ratio of K10_4, K10_5, KA5_2 by mass spec
def t10_ms():
    i = 7
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    span = [23.6, 24.0]
    label="K10_4,K10_5,KA5_2"

    span_ms = ms.extract_ms(span[0], span[1])
    fig = ca.plotting.signal(span_ms, normalize=1)


    lib = ca.library.Library.from_JSON(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.json")
    ms_raw = lib.find_by_name("4-decanone").get_attribute("MassSpectrum").value
    K10_4 = ca.mass_spec.MSSignal(x=ms_raw[:,0], y=ms_raw[:,1], name="K10_4")
    fig = ca.plotting.signal(K10_4, fig=fig, normalize=1)

    ms_raw = lib.find_by_name("5-decanone").get_attribute("MassSpectrum").value
    K10_5 = ca.mass_spec.MSSignal(x=ms_raw[:,0], y=ms_raw[:,1], name="K10_5")
    fig = ca.plotting.signal(K10_5, fig=fig, normalize=1)

    ms_raw = np.round(np.array(KA5_2_ms_raw, dtype=float))
    KA5_2 = ca.mass_spec.MSSignal(x=ms_raw[:,0], y=ms_raw[:,1], name="KA5_2")
    fig = ca.plotting.signal(KA5_2, fig=fig, normalize=1)

    sigs = ca.mass_spec.unify_mz((span_ms, K10_4, K10_5, KA5_2))
    X = np.column_stack([sigs[1].y_normalized_by_max(), sigs[2].y_normalized_by_max(), sigs[3].y_normalized_by_max()])   # shape (N, 3)


    from scipy.optimize import nnls
    w, rnorm = nnls(X, sigs[0].y_normalized_by_max())

    print(w,w/np.sum(w), rnorm, rnorm / np.linalg.norm(sigs[0].y_normalized_by_max()))
    #[0.30247269 0.16452585 0.94635134] [0.21401119 0.11640844 0.66958037] 0.2702942118989266 0.1594587450516134
    # fig.show("browser")


    fit = X @ w
    fig2 = ca.plotting.signal(span_ms, normalize=1)
    fig2.add_bar(x=sigs[0].x, y=fit, name="fit")
    fig2.show("browser")



def plot_series():
    fig = go.Figure()
    colors = plotly_utils.get_colors_from_scale(len(raw_data), "Viridis")
    colors.reverse()
    for i, data in enumerate(raw_data):
        x, y = data[1].x, data[1].y_normalized_by_max([28, 30])
        y = y  # + 0.02*i
        fig.add_scatter(x=x, y=y, mode="lines", line_color=colors[i], name=f"{times[i]} min")

    # full
    # fig.layout.xaxis.range = 5.5, None
    # fig.layout.yaxis.range = -0.5, 7

    # zoom 1
    # fig.layout.xaxis.range = 5.5, 25
    # fig.layout.yaxis.range = -0.02, 1

    # zoom 2
    # fig.layout.xaxis.range = 20, 45
    # fig.layout.yaxis.range = -0.02, 1

    # diketone
    fig.layout.xaxis.range = 38, 45
    fig.layout.yaxis.range = -0.01, 0.1

    fig.layout.xaxis.title.text = "retention time (min)"
    fig.layout.yaxis.title.text = "FID intensity"
    fig.layout.margin.t = 20
    # fig.show('browser')
    fig.write_image("DJW-11-115-V1_diketone.png")  # _zoom1


def main():
    data = [t0(), t1(), t2(), t3(), t4(), t5(), t6(), t7(), t8(), t9(), t10(), t11()]  # t12(), t13(), t14()
    df = pl.DataFrame(data)
    df.with_columns(pl.col(col).fill_null(0) for col in df.columns)
    print(df)
    df.write_csv(data_label + "_areas.csv")

    # print(t6())


if __name__ == "__main__":
    # plot_series()
    # main()
    t10_ms()
