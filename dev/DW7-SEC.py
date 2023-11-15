
import numpy as np

import chem_analysis as ca


def main2():
    def cal_func(x: np.ndarray) -> np.ndarray:
        return 10**(-0.6035623045394646*x + 10.70478909408625)
    cal = ca.sec.ConventionalCalibration(cal_func, time_bounds=[8.377, 13.1])
    data = ca.sec.SECSignalArray.from_file(
        r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_SEC.feather",
        calibration=cal
    )

    stats = ""
    for i in range(data.time.size): # data.x.size
        signal = data.get_signal(i)
        result_peaks_max = ca.analysis.peak_picking.scipy_find_peaks(signal, height=2, distance=500, prominence=1)
        if len(result_peaks_max.indexes) == 0:
            print(i, "no peaks")
            stats += f"\n{i} No peaks"
            continue
        result_peaks = ca.analysis.boundary_detection.rolling_ball(result_peaks_max, n=15)
        if stats == "":
            stats += "signal," + ",".join(result_peaks.stats_table().headers)

        stats += f"\n{i}," + result_peaks.stats_table().to_csv(with_headers=False)
        print(i, len(result_peaks.peaks))

    with open("table.csv", "w", encoding="UTF-8") as f:
        f.write(stats)


def main():
    def cal_func(x: np.ndarray) -> np.ndarray:
        return 10**(-0.6035623045394646*x + 10.70478909408625)
    cal = ca.sec.ConventionalCalibration(cal_func, time_bounds=[8.377, 13.1])
    data = ca.sec.SECSignalArray.from_file(
        r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_SEC.feather",
        calibration=cal
    )

    # single
    signal = data.get_signal(12)
    signal.processor.add(
        ca.processing.baseline_correction.Polynomial(
            degree=2,
            mask=ca.processing.masks.DataMaskDistanceRemove(amount=0.2)
        )
    )
    result_peaks_max = ca.analysis.peak_picking.scipy_find_peaks(signal, height=2, distance=500, prominence=1)
    result_peaks = ca.analysis.boundary_detection.rolling_ball(result_peaks_max, n=20, n_points_with_pos_slope=5)
    print(result_peaks.stats_table().to_str(sig_figs=4))

    # with open("table.csv", "w", encoding="UTF-8") as f:
    #     f.write(result_peaks.stats_table().to_csv(with_headers=False))

    fig = ca.sec.plot_sec_calibration(signal.calibration)
    fig = ca.sec.plot_sec_peaks(result_peaks, fig=fig)
    fig = ca.sec.plot_sec_signal(signal, fig=fig)
    fig.write_html("temp.html", auto_open=True)
    # print("hi")

def main3():
    from pybaselines import Baseline, utils

    def cal_func(x: np.ndarray) -> np.ndarray:
        return 10**(-0.6035623045394646*x + 10.70478909408625)
    cal = ca.sec.ConventionalCalibration(cal_func, time_bounds=[8.377, 13.1])
    data = ca.sec.SECSignalArray.from_file(
        r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_SEC.feather",
        calibration=cal
    )

    # single
    signal = data.get_signal(12)

    mask= ca.processing.masks.DataMaskDistanceRemove(amount=0.5, keep_ends=True)
    x, y = mask.get_data(signal.x, signal.y)

    baseline_fitter = Baseline(x_data=x)

    bkg_1 = baseline_fitter.modpoly(y, poly_order=3)[0]
    bkg_2 = baseline_fitter.asls(y, lam=1e7, p=0.02)[0]

    import plotly.graph_objs as go
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=signal.x, y=signal.y))
    fig.add_trace(go.Scatter(x=x, y=bkg_2))
    fig.write_html("temp.html", auto_open=True)


if __name__ == "__main__":
    # main()
    # main2()
    main3()
