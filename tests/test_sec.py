import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca


def calibration_func(retention_time: np.ndarray) -> np.ndarray:
    return 10 ** (-0.3167*retention_time+11.294)  # mw_i


def calibration_func_inverse(mw_i: np.ndarray) -> np.ndarray:
    return (11.294 - np.log10(mw_i))/0.3167  # retention time


def log_normal_wi(mw_i: np.ndarray, Mn: float = 60_000, D: float = 1.5) -> np.ndarray:
    weight_fraction = 1 / (Mn * np.sqrt(2 * np.pi * np.log(D))) * \
                      np.exp(-1 * (np.log(mw_i / Mn) + np.log(D) / 2) ** 2 / (2 * np.log(D)))
    return weight_fraction


def compute_mn_d(signal: ca.sec.SECSignal):
    # signal.processor.add(ca.processing.baseline_correction.Polynomial(degree=3))
    peaks = ca.analysis.peak_picking.max_find_peaks(signal, weights=ca.processing.weigths.Spans((13, 30), invert=True))
    peaks = ca.analysis.boundary_detection.rolling_ball(peaks, n=10, min_height=0.0001, n_points_with_pos_slope=1)
    print(peaks.stats_table().write_csv_str())

    fig = ca.plot.calibration(signal.calibration)
    fig = ca.plot.peaks(peaks, fig=fig)
    fig = ca.plot.signal(signal, fig=fig)

    fig.layout.yaxis.range = (-1, peaks.max_height)
    fig.show()


def main():
    Mn = 60_000
    D = 1.6
    n = 10000
    mw_i = np.linspace(100, Mn*10, n)
    w_i = log_normal_wi(mw_i=mw_i, Mn=Mn, D=D)
    retention_time = calibration_func_inverse(mw_i)
    x_i = w_i * Mn / mw_i

    # fig = go.Figure()
    # fig.add_trace(go.Scatter(x=retention_time, y=w_i))
    # fig.show()

    calibration = ca.sec.ConventionalCalibration(calibration_func, mw_bounds=(104, 1_000_000))
    signal = ca.sec.SECSignal(x_raw=retention_time, data_raw=w_i, calibration=calibration)
    compute_mn_d(signal)


if __name__ == "__main__":
    main()
