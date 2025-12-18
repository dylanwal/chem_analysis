import pathlib

import numpy as np
import polars as pl
import plotly.graph_objs as go

import chem_analysis as ca

import development.plotly_utils as plotly_utils

folder_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\11\11_115\GCMS")
data_label = "DJW-11-115-V1-t"
times = np.array([0, 30, 60, 90, 120, 150, 180, 240, 300, 360, 480, 600, 1440, 2160, 2880])


def load_single_file(data_path: pathlib.Path, label: str) -> tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]:
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path/ f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    return ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)

labels = [data_label + str(i) for i in times]
raw_data: list[tuple[ca.gc_lc.GCMSSignal, ca.gc_lc.GCSignal]] = [load_single_file(folder_path, label) for label in labels if "600" not in label]

def load_600() -> tuple[ca.gc_lc.GCSignal, ca.gc_lc.GCSignal]:
    ms_data = np.loadtxt(r"C:\Users\nicep\Desktop\research_wis\data\11\11_115\GCMS\DJW-11-115-V1-t600_ms.csv", delimiter=",", skiprows=2)
    ms = ca.gc_lc.GCSignal(x=ms_data[:, 1], y=ms_data[:, 2])

    fid_data = np.loadtxt(r"C:\Users\nicep\Desktop\research_wis\data\11\11_115\GCMS\DJW-11-115-V1-t600_fid.csv", delimiter=",", skiprows=2)
    fid = ca.gc_lc.GCSignal(x=fid_data[:, 1], y=fid_data[:, 2])

    return ms, fid

raw_data.insert(11, load_600())


def plot_ms_fid(i: int, ms: ca.gc_lc.GCMSSignal, fid: ca.gc_lc.GCSignal, peaks: list[ca.a.PeakData]):
    fig_ms = ca.plotting.signal(ms)

    fig_fid = ca.plotting.signal(fid)
    fig_fid = ca.plotting.add_peaks(peaks, fig=fig_fid, colors="multi")
    fig_fid.layout.yaxis.range = [-50_000, 4_000_000]

    plotly_utils.merge_html_figs([fig_ms, fig_fid], filename=f"plots_{times[i]}.html", auto_open=True)


def t0():
    i=0
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [28.96, 29.64], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True)
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t1():
    i=1
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.8], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True)
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t2():
    i=2
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.53], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [8.6, 9.25], label="decane", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [9.8, 10.15], label="CA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [14.9, 15.3], label="CA6", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [15.3, 15.75], label="HA2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [20.8, 21.3], label="CA7", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [23.6, 24.0], label="K10_4,K10_5,KA5_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.25, 24.4], label="A10_4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.4, 24.7], label="HAA3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [24.8, 25.0], label="?HAA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t3():
    i=3
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t4():
    i=4
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t5():
    i=5
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t6():
    i=6
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t7():
    i=7
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t8():
    i=8
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t9():
    i=9
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t10():
    i=10
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t11():
    i=11
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data


def t12():
    i=12
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t13():
    i=13
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def t14():
    i=14
    data = dict()
    data["time"] = times[i]
    ms, fid = raw_data[i]
    peaks =  [
        ca.a.peaks.peak_from_span(fid,  [6.2, 6.48], label="CA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [8.30, 8.48], label="K7_2", to_zero=True),
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
        ca.a.peaks.peak_from_span(fid, [25.0, 25.15], label="K10_3,", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [25.5, 25.84], label="K10_2,A10_3", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [26.6, 26.85], label="A2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [27.06, 27.33], label="CA8", to_zero=True),
        ca.a.peaks.peak_from_span(fid,  [28.94, 29.32], label="TCB", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [29.32, 29.51], label="KA6_3", to_zero=True),
        # ca.a.peaks.peak_from_span(fid, [29.68, 29.9], label="?", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [30.3, 30.56], label="KA6_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32, 32.3], label="DA4", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.6, 32.8], label="DA4_double", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [32.8, 33.0], label="AA10_2", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [37.6, 37.9], label="DA5", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [38.9, 39.2], label="CA10", to_zero=True),
        ca.a.peaks.peak_from_span(fid, [43.6, 43.77], label="DA5", to_zero=True),  #mix
    ]

    plot_ms_fid(i, ms, fid, peaks)
    data.update({peak.label: peak.properties.area for peak in peaks})
    return data

def plot_series():
    fig = go.Figure()
    colors = plotly_utils.get_colors_from_scale( len(raw_data),"Viridis")
    for i, data in enumerate(raw_data):
        x, y = data[1].x, data[1].y_normalized_by_max([28, 30])
        y = y #+ 0.02*i
        fig.add_scatter(x=x, y=y, mode="lines", line_color=colors[i], name=f"{times[i]} min")

    fig.layout.xaxis.range = 5.5, None
    fig.layout.xaxis.title.text = "time (min)"
    fig.layout.yaxis.title.text = "FID intensity"
    fig.layout.yaxis.range = -0.5, 7
    # fig.layout.yaxis.range = -0.02, 0.5
    fig.layout.margin.t = 20
    fig.show('browser')


def main():
    # data = [t0(), t1(), t2(), t3(), t4(), t5(), t6(), t7(), t8(), t9(), t10(), t11(), t12(), t13(), t14()]
    # df = pl.DataFrame(data)
    # df.with_columns(pl.col(col).fill_null(0) for col in df.columns)
    # print(df)
    # df.write_csv(data_label + "_areas.csv")

    print(t6())


if __name__ == "__main__":
    # plot_series()
    main()