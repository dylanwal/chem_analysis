import pathlib

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
    ms_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(ms, fig=ms_fig)
    ca.plotting.peaks(ms_peaks, fig=ms_fig)
    ms_fig.layout.title = data_path.stem
    ms_fig.layout.title = "MS"

    ms_ms = ca.analysis.ms_extraction.ms_extract_peak(ms_peaks[0])
    ms_ms.processor.add(ca.processing.edit.CutSpans([None, 5]))
    ms_ms_fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ms_ms_fig.layout.title = data_path.stem
    ms_ms_fig.layout.xaxis.range = (10, 60)
    ca.plotting.signal(ms_ms, fig=ms_ms_fig)

    print(ms_ms.to_numpy(reduce=True, cutoff=0.01))

    print("finished analyzing:", data_path.name)
    return [ms_fig, ms_ms_fig]


def main():
    folder_path = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\gas")
    folders = [item for item in folder_path.iterdir() if item.is_dir()]
    figs = map(process_single, folders)

    figs = [item for sublist in figs for item in sublist]
    ca.plotting.PlotlyConfig.merge_figures(figs, filename=folder_path / "ms", auto_open=True)


if __name__ == '__main__':
    main()
