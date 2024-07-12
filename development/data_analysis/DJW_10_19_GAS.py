import pathlib
import re

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca


def process_single(data_path: pathlib.Path):
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_folder(data_path)
    # ms
    ms.processor.add(
        ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4, save_result=True)
    )
    peak_locations = ca.analysis.peak_picking.find_peaks_scipy(ms,
                                                               scipy_kwargs={"prominence": 10000}
                                                               )
    ms_peaks = ca.analysis.integration.rolling_ball(peak_locations, n=5, min_height=0.002,
                                                    n_points_with_pos_slope=2)

    # plotting peak results
    gc_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(ms, fig=gc_fig)
    ca.plotting.peaks(ms_peaks, fig=gc_fig)
    gc_fig.layout.title = data_path.stem
    gc_fig.layout.title = "MS"

    ms_ms = ca.analysis.ms_extraction.ms_extract_peak(ms_peaks[0])
    print("finished analyzing:", data_path.name)
    return ms, ms_ms, gc_fig


def get_folders(data_path: pathlib.Path, pattern: str) -> list[pathlib.Path]:
    pattern_compile = re.compile(pattern.replace("*", "([0-9]+)"))

    def sort_func(x: pathlib.Path) -> int:
        return int(pattern_compile.match(x.name).groups()[0])

    specific_folders = list(data_path.glob(pattern))
    print(len(specific_folders), "files found for analysis")

    specific_folders.sort(key=sort_func)

    if len(specific_folders) == 0:
        raise RuntimeError("No folders found for analysis")

    return [data_path / file for file in specific_folders]


def intensity_to_mmol(CO2: int| float, Ar: int | float) -> int | float:
    area_ratio = CO2/Ar
    molar_ratio = 44.6*area_ratio**0.767

    return molar_ratio*0.33  # mmol/min


def process_timeseries(data_path: pathlib.Path | str, pattern: str):
    if isinstance(data_path, str):
        data_path = pathlib.Path(data_path)
    folders = get_folders(data_path, pattern)

    figs = []
    ms = []
    ms_ms = []
    times = np.empty(len(folders))
    for i, folder in enumerate(folders):
        ms_, ms_ms_, fig = process_single(folder)
        ms.append(ms_)
        ms_ms.append(ms_ms_)
        times[i] = ms_.parameters.time_start_run.timestamp()
        if i == int(len(folders)/2):
            figs.append(fig)
            ms_ms_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
            ms_ms_fig.layout.title = data_path.stem
            ms_ms_fig.layout.xaxis.range = (10, 60)
            ca.plotting.signal(ms_ms_, fig=ms_ms_fig)
            figs.append(ms_ms_fig)

    # plot timeseries
    times -= times[0]
    times /= 60
    times -= 30
    y = np.empty_like(times)
    for i in range(len(times)):
        CO2 = ms_ms[i].get_intensity(44)
        Ar = ms_ms[i].get_intensity(40)
        y[i] = intensity_to_mmol(CO2, Ar)
    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    fig.add_scatter(x=times, y=y)
    fig.layout.xaxis.title = "<b>time (min) </b>"
    fig.layout.yaxis.title = "<b>mmol/min </b>"
    fig.update_layout(
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)'
)
    figs.append(fig)
    ca.plotting.PlotlyConfig.merge_figures(figs, filename=data_path / "GAS_analysis", auto_open=True)
    print("total CO2 mmol:", np.trapz(x=times, y=y))


def main():
    data_path = r"C:\Users\nicep\Desktop\research_wis\data\10\10_19\GC_MS"
    pattern = "DJW-19-GAS-*.D"
    process_timeseries(data_path, pattern)


if __name__ == '__main__':
    main()
