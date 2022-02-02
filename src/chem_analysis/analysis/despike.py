import numpy as np
from numba import njit, config

config.DISABLE_JIT = True


@njit
def detect_outliers(data: np.ndarray, m: float = 2):
    dist_from_median = np.abs(data - np.median(data))
    median_deviation = np.median(dist_from_median)
    if median_deviation != 0:
        scale_distances_from_median = dist_from_median/median_deviation
        return scale_distances_from_median > m  # True is an outlier

    return np.zeros_like(data)  # no outliers


# @njit
def window_calc(data: np.ndarray, pos: int, m: float = 2) -> float:
    outlier = detect_outliers(data, m)
    if outlier[pos]:
        return np.median(data[np.invert(outlier)])
    else:
        return data[pos]


# @njit
def despike(y: np.ndarray, window: int = 20, m: float = 2) -> np.ndarray:
    out = np.empty_like(y)

    if window % 2 != 0:
        window += 1

    span = int(window/2)
    for i in range(len(y)):
        if i < span:  # left edge
            out[i] = window_calc(y[:window], i, m)
        elif i > len(y) - span:  # right edge
            out[i] = window_calc(y[window:], i - (len(y)-window), m)
        else:  # middle
            out[i] = window_calc(y[i-span:i+span], span, m)

    return out


def test():
    import plotly.graph_objs as go

    n = 1000
    x = np.linspace(0, n-1, n)
    y = np.ones(n) + np.random.random(n)
    y[100] = 5
    y[101] = 3

    y[500] = 2
    y[501] = 1.7
    y[502] = 1.5
    y[503] = 1.2

    y[700] = 0.5
    y[701] = 0
    y[702] = -0.5
    y[703] = -1

    y_new = despike(y)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode="lines"))
    fig.add_trace(go.Scatter(x=x, y=y_new, mode="lines"))
    fig.write_html("temp.html", auto_open=True)


if __name__ == "__main__":
    test()
