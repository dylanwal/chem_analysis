
import numpy as np

import chem_analysis as ca


def main():
    def cal_func(x: np.ndarray) -> np.ndarray:
        return 10**(-0.6035623045394646*x + 10.70478909408625)
    cal = ca.sec.ConventionalCalibration(cal_func, time_bounds=[8.377, 13.1])

    data = ca.sec.SECSignalArray.from_file(
        r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_SEC.feather",
        calibration=cal
    )

    signal = data.get_signal(0)
    result_peaks = ca.analysis.peak_picking.scipy_find_peaks(signal, height=20, distance=300)
    result_peaks.indexes = [result_peaks.indexes[0]]
    result_peaks = ca.analysis.boundry_detection.rolling_ball(result_peaks)
    print(result_peaks.stats_table().to_str(sig_figs=4))
    fig = ca.sec.plot_sec_signal(signal)
    fig = ca.sec.plot_sec_calibration(signal.calibration, fig=fig)
    fig = ca.sec.plot_sec_peaks(result_peaks, fig=fig)
    fig.write_html("temp.html", auto_open=True)
    print("hi")


if __name__ == "__main__":
    main()
