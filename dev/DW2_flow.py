from typing import Callable
from collections import OrderedDict
import functools

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import plotly.graph_objs as go

from chem_analysis.utils.general_math import get_slice
from chem_analysis.utils.feather_format import feather_to_numpy, unpack_time_series


def distribution_normal(x, amp=1, mean=0, sigma=1):
    return amp / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))


# from scipy.stats import cauchy
def distribution_cauchy(x, amp=1, gamma=1, mean=0):
    return amp / (np.pi * gamma * (1 + ((x - mean) / 2) ** 2))


from scipy.special import voigt_profile
def distribution_volt(x, amp, sigma, gamma, mean):
    """
    gamma = 0 normal
    sigma = 0 cauchy distribution
    """
    x = x - mean
    return amp * voigt_profile(x, sigma=sigma, gamma=gamma)


def multi_distribution(func: Callable, x: np.ndarray, n: int, *args) -> np.ndarray:
    if (len(args) % n) != 0:
        raise ValueError("issue with args")

    y = np.zeros_like(x, dtype=x.dtype)
    arg_n = int(len(args)/n)
    for i in range(n):
        y += func(x, args[arg_n*i:arg_n*(i+1)])

    return y


def multi_distribution_integration(func: Callable, x: np.ndarray, n: int, args) -> np.ndarray:
    if (len(args) % n) != 0:
        raise ValueError("issue")

    area = np.zeros_like(n, dtype=x.dtype)
    arg_n = int(len(args)/n)
    for i in range(n):
        y = func(x, args[arg_n*i:arg_n*(i+1)])
        area[i] = np.trapz(x=x, y=y)

    return area


def conv_from_integration(data, range_1: list, range_2: list):
    slice_1 = get_slice(data.x, range_1[0], range_1[1])
    area_1 = np.trapz(x=data.x[slice_1], y=data.data[:, slice_1])

    slice_2 = get_slice(data.x, range_2[0], range_2[1])
    area_2 = np.trapz(x=data.x[slice_2], y=data.data[:, slice_2])

    return area_2 / (area_1 + area_2)


def get_peaks(data, slice_, n: int = 10):
    prob_peaks = np.zeros((n, 2), dtype=data.data.dtype)
    for i in range(n):
        ii = np.random.randint(0, len(data.time))
        y = data.data[ii, slice_]
        peaks, prop = find_peaks(y, prominence=1, width=0.1)
        peaks[i, :] = peaks

    return np.mean(prob_peaks, axis=0)


def conv_from_normal(data):
    range_ = [3.4, 3.9]
    p0 = OrderedDict(
        (
            ("amp1", 1),
            ("mean1", 3.53),
            ("sigma1", 0.01),
            ("amp2", 1),
            ("mean2", 3.65),
            ("sigma2", 0.01)
        )
    )
    slice_ = get_slice(data.x, start=range_[0], end=range_[1])
    x = data.x[slice_]

    areas = np.empty((data.time.size, 2), dtype=np.float64)
    issues = []
    good = []
    for i, row in enumerate(data.data):
        try:
            result = curve_fit(
                functools.partial(multi_distribution, func=distribution_normal, n=2),
                x,
                data.data[i, slice_],
                p0=list(p0.values())
            )
            params, covariance = result
            areas[i] = multi_distribution_integration(func=distribution_normal, n=2, x=data.x, args=params)
        except Exception as e:
            issues.append(i)
            print(e)
            continue

        good.append(i)

    print("issues: ", len(issues), "| total:", len(data.time))
    areas = areas[good]
    times_ = data.time[good]

    return np.column_stack((times_, areas[:, 1] / (areas[:, 0] + areas[:, 1])))


class NMRData:
    def __init__(self, ppm, times, data):
        self.time = times
        self.time_zeroed = self.time - self.time[0]
        self.pmm = ppm
        self.data = data

    @property
    def x(self):
        return self.pmm


def conv_analysis(data: NMRData):
    conv_integration = conv_from_integration(data, range_1=[3.65, 3.8], range_2=[3.488, 3.65])
    print("integration:", conv_integration[-1])
    conv_normal = conv_from_normal(data)
    # conv_cauchy = conv_from_fit(data, range_=range_, func=distribution_cauchy)
    # conv_voigt = conv_from_fit(data, range_=range_, func=distribution_volt)


    fig = go.Figure()
    fig.add_trace(go.Scatter(x=data.x, y=conv_integration, name="conv_int"))
    # fig.add_trace(go.Scatter(x=data.x, y=conv_normal, name="conv_normal"))
    # fig.add_trace(go.Scatter(x=data.x, y=conv_cauchy, name="conv_cauchy"))
    # fig.add_trace(go.Scatter(x=data.x, y=conv_voigt, name="conv_cauchy"))
    fig.update_layout(autosize=False, width=800, height=600, font=dict(family="Arial", size=18, color="black"),
                      plot_bgcolor="white", showlegend=True)
    fig.update_xaxes(title="<b>rxn time (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray")
    fig.update_yaxes(title="<b>conversion</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray", range=[0, 1])
    fig.show()


def plot(data):
    n = [0, 2, 5, 10, 20, 30, 200, -1]
    fig = go.Figure()
    for i in n:
        fig.add_trace(go.Scatter(x=data.x, y=data.data[i]))

    fig.show()


def clean_data(data):
    remove_spectra = []
    for i in range(len(data.time)):
        if np.max(data.data[i]) < 10:
            remove_spectra.append(i)

    data.time = np.delete(data.time, remove_spectra)
    data.time_zeroed = np.delete(data.time_zeroed, remove_spectra)
    data.data = np.delete(data.data, remove_spectra, axis=0)


def main():
    path = r"C:\Users\nicep\Desktop\DW2_flow_rate_NMR.feather"
    data = feather_to_numpy(path)
    data = NMRData(*unpack_time_series(data))
    clean_data(data)

    # plot(data)
    conv_analysis(data)


if __name__ == "__main__":
    main()
