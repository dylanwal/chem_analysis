
import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca


def get_time_series(folder_path: str, pattern: str) -> ca.gc_lc.GCMSSignal2D:
    import glob

    files = glob.glob(folder_path + "/" + pattern)
    files.sort(key=lambda x: int(x[66:-9]))
    times = np.array([int(f[66:-9]) for f in files])
    data = [ca.gc_lc.GCParser.from_Agilent_D_folder(f)[0] for f in files]

    return ca.gc_lc.GCMSSignal2D.from_signals(data, times)


def first_load(folder_: str, pattern: str):
    data = get_time_series(folder_, pattern)
    print(data)
    data.to_npy(folder_ + r"\data")


def main(folder_: str):
    lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.JSON"
    chemistry_lib = ca.mass_spec.GCLibrary.from_JSON(lib_path)

    data = ca.gc_lc.GCMSSignal2D.from_file(folder_ + r"\data.npz")
    data = data.get_signal(3)
    data.processor.add(ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4))

    # signal_ = data.get_signal(3)
    # fig = go.Figure()
    # fig.add_trace(go.Scatter(x=signal_.x, y=signal_.y))
    # fig.show()

    peaks = ca.analysis.peak_picking.find_peaks_scipy(data, scipy_kwargs={"height": 5000, "width": 0.1})
    peak_result = ca.analysis.boundary_detection.rolling_ball(peaks, n=10, min_height=0.05, n_points_with_pos_slope=1)

    picking_lib = chemistry_lib.to_picking_library()
    peaks_results_identified = ca.analysis.peak_picking.library_search.find_peaks_retention_time_library(peak_result, picking_lib)
    print(peaks_results_identified)


if __name__ == "__main__":
    folder = r"C:\Users\nicep\Desktop\research_wis\data\11\11_23\gc_ms"
    # first_load(folder, "DJW-11-23-*min-TMS.D")
    main(folder)
