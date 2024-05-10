import glob
import pathlib
from datetime import datetime
from collections import OrderedDict

import numpy as np
import plotly.graph_objs as go

import chem_analysis.utils.math
from chem_analysis.gc_lc.parsers.agilent_folder import parse_D_folder

# from chem_analysis.processing.baseline import bc_polynomial

template = go.layout.Template()
template.layout.font = dict(family="Arial", size=18, color="black")
template.layout.plot_bgcolor = "white"
template.layout.width, template.layout.height = 1200, 600
template.layout.xaxis.tickprefix = "<b>"
template.layout.xaxis.ticksuffix = "<b>"
template.layout.xaxis.showline = True
template.layout.xaxis.linewidth = 5
template.layout.xaxis.linecolor = "black"
template.layout.xaxis.ticks = "outside"
template.layout.xaxis.tickwidth = 4
template.layout.xaxis.showgrid = False
template.layout.xaxis.mirror = True
template.layout.yaxis.tickprefix = "<b>"
template.layout.yaxis.ticksuffix = "<b>"
template.layout.yaxis.showline = True
template.layout.yaxis.linewidth = 5
template.layout.yaxis.linecolor = "black"
template.layout.yaxis.ticks = "outside"
template.layout.yaxis.tickwidth = 4
template.layout.yaxis.showgrid = False
template.layout.yaxis.mirror = True


def gas_ratio_single(data: tuple):
    ini_dict, ms_dict, fid_dict = data

    ms_data = ms_dict['data']
    ms_time = ms_dict['time']
    ms_sum = np.sum(ms_data, axis=1)
    co2 = ms_data[:, 44]
    o2 = ms_data[:, 32]
    x_range = chem_analysis.utils.math.get_slice(ms_time, 1.7568, 2.130)
    O2_area = np.trapz(x=ms_time[x_range], y=o2[x_range])
    CO2_area = np.trapz(x=ms_time[x_range], y=co2[x_range])
    x_range = chem_analysis.utils.math.get_slice(ms_time, 5.4, 5.48)
    decane_area = np.trapz(x=ms_time[x_range], y=ms_sum[x_range])

    time_: datetime = ms_dict['date_time']
    return time_.timestamp(), O2_area, CO2_area, decane_area


def gas_ratio(data):
    O2 = np.empty(len(data))
    CO2 = np.empty_like(O2)
    decane = np.empty_like(O2)
    time = np.empty_like(O2)
    for i, f in enumerate(data):
        time[i], O2[i], CO2[i], decane[i] = gas_ratio_single(f)
        # print(i, time[i], O2[i], CO2[i], decane[i])

    return time, O2, CO2, decane


def main_gas():
    folder = r"C:\Users\nicep\Desktop\research_wis\data\11\11_23\gc_ms"
    files = glob.glob(folder + "/" + "DJW-11-23-GAS-*.D")
    files.sort(key=lambda x: int(x[70:-2]))
    data = [parse_D_folder(pathlib.Path(f)) for f in files]
    time, O2, CO2, decane = gas_ratio(data)

    time_start = datetime.fromisoformat('2024-05-01T10:31:00').timestamp()
    time_o2 = datetime.fromisoformat('2024-05-01T09:58:00').timestamp() - time_start
    time_normal = time - time_start

    fig = go.Figure(layout={"template": template})
    fig.add_trace(go.Scatter(x=time_normal[1:] / 60, y=O2[1:], mode='lines+markers', name="O2"))
    fig.add_trace(go.Scatter(x=time_normal[1:] / 60, y=CO2[1:], mode='lines+markers', name="CO2"))
    fig.add_trace(go.Scatter(x=time_normal[1:] / 60, y=decane[1:], mode='lines+markers', name="decane"))
    fig.add_trace(go.Scatter(x=[time_o2 / 60, time_o2 / 60], y=[0, np.max(O2)], mode='lines', name="O2 start"))
    fig.add_trace(go.Scatter(x=[0, 0], y=[0, np.max(O2)], mode='lines', name="start heat"))
    fig.layout.xaxis.title = "<b>time (min)</b>"
    fig.layout.yaxis.title = "<b>ion count</b>"
    fig.layout.yaxis.range = [0, None]
    fig.show()

    fig = go.Figure(layout={"template": template})
    fig.add_trace(go.Scatter(x=data[15][1]['time'], y=np.sum(data[15][1]['data'], axis=1)))
    fig.show()


#######################################################################################################################

chems = OrderedDict(
    CA3=(3.887, 4.059),
    CA4=(5.95, 6.244),
    C10=(9, 9.24),
    CA5=(9.96, 10.219),
    CA6=(15.312, 15.542),
    HAA2=(17.81, 18.10),
    CA7=(21.36, 21.65),
    K10_4=(24.13, 24.338),
    A10_5=(24.338, 24.553),
    A10_4=(24.8467, 25.118),
    K10_3=(25.39, 25.61),
    K10_2=(25.96, 26.18),
    A10_3=(26.18, 26.4),
    A10_2=(27.117, 27.382),
    CA8=(27.61, 27.86),
    TCB=(29.5, 30.0),
    DA4=(32.6, 33),
    DA5=(38.2, 38.7)
)


def liquid_single(data):
    ini_dict, ms_dict, fid_dict = data

    ms_data = ms_dict['data']
    ms_time = ms_dict['time']
    ms_sum = np.sum(ms_data, axis=1)
    baseline = baseline_section_std(ms_sum, sections=100, window=15, number_of_deviations=4)
    ms_sum = ms_sum - baseline

    areas = np.empty(len(chems))
    for i, chem in enumerate(chems):
        x_range = chem_analysis.utils.math.get_slice(ms_time, *chems[chem])
        areas[i] = np.trapz(x=ms_time[x_range], y=ms_sum[x_range])

    return areas


def liquid(data):
    results = np.empty((len(data), len(chems)))
    for i, f in enumerate(data):
        results[i, :] = liquid_single(f)

    return results


def main_liq():
    folder = r"C:\Users\nicep\Desktop\research_wis\data\11\11_23\gc_ms"
    files = glob.glob(folder + "/" + "DJW-11-23-*min-TMS.D")
    files.sort(key=lambda x: int(x[66:-9]))
    # files = [files[2], files[9]]
    times = np.array([int(f[66:-9]) for f in files])
    times = times - times[0]
    data = [parse_D_folder(pathlib.Path(f)) for f in files]
    results = liquid(data)

    fig = go.Figure(layout={"template": template})
    for i, chem in enumerate(chems):
        fig.add_trace(go.Scatter(x=times, y=results[:, i], name=chem))
    fig.add_trace(go.Scatter(x=times, y=np.sum(results, axis=1), name="mass balance"))
    fig.layout.xaxis.type = "log"
    fig.show()


def main_liq_baseline():
    folder = r"C:\Users\nicep\Desktop\research_wis\data\11\11_23\gc_ms"
    files = glob.glob(folder + "/" + "DJW-11-23-*min-TMS.D")
    files.sort(key=lambda x: int(x[66:-9]))
    files = [files[2], files[9]]
    data = [parse_D_folder(pathlib.Path(f)) for f in files]
    # liquid(data)

    # fig = go.Figure(layout={"template": template})
    # fig.add_trace(go.Scatter(x=data[1][1]['time'], y=np.sum(data[1][1]['data'], axis=1)))
    # fig.show()

    x = data[0][1]['time']
    y = np.sum(data[0][1]['data'], axis=1)
    baseline = baseline_section_std(y, sections=100, window=15, number_of_deviations=4)
    fig = go.Figure(layout={"template": template})
    fig.add_trace(go.Scatter(x=x, y=y))
    fig.add_trace(go.Scatter(x=x, y=baseline))
    fig.add_trace(go.Scatter(x=x, y=y - baseline))
    fig.layout.yaxis.range = -15000, 15000
    fig.show()


##########################################################################################################

def baseline_section_std(y,
                         window: int = 3,
                         sections: int = 32,
                         number_of_deviations: int | float = 2,
                         smooth_window: int = 15,
                         ):
    mask = sectioned_std(y, window, sections, number_of_deviations, smooth_window)
    x = np.arange(len(y))
    mask[0], mask[-1] = True, True
    x_mask = x[mask]
    y_mask = y[mask]
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(np.interp(x, x_mask, y_mask), 100)


def sectioned_std(y,
                  window: int = 3,
                  sections: int = 32,
                  number_of_deviations: int | float = 2,
                  smooth_window: int = 20,
                  ):
    """
    This algorithm assumes the y data contains at least one region with no signals.

    https://doi.org/10.1006/jmre.2000.2121

    Parameters
    ----------
    y:
        data
    window:
        number of points used to compute local variation
    sections:
        number of sections used to get minimum standard deviation || noise value
    number_of_deviations:
        the number of deviations from the noise value
        little effect, typically 2 to 4
    smooth_window:
        little effect
        smoothing window size

    Returns
    -------
    mask:
        1 where the baseline is
        0 where peaks are
    """
    from chem_analysis.utils.pad_edges import pad_edges_polynomial
    from numpy.lib.stride_tricks import sliding_window_view

    window = max(window, 1)
    # compute noise level by breaking the data into sections and find section with min sigma
    slices = divide_array(y, sections)
    min_sigma = np.inf
    for slice_ in slices:
        min_sigma = min(min_sigma, np.std(y[slice_]))

    # smooth specta with convolution with rectangular function
    # half_window = int(smooth_window/2)
    # padded_y = pad_edges_polynomial(y, pad_amount=half_window)
    # sliding_window = sliding_window_view(padded_y, 2*half_window+1)
    # smoothed_y = np.sum(sliding_window, axis=1) / (2*half_window+1)

    from scipy.ndimage import gaussian_filter
    smoothed_y = gaussian_filter(y, 5)

    # evaluate if point is outside min_sigma
    # (max(y_i) - min(y_i)) < n*sigma
    half_window = int(window / 2)
    padded_smoothed_y = pad_edges_polynomial(smoothed_y, pad_amount=half_window)
    sliding_window = sliding_window_view(padded_smoothed_y, 2 * half_window + 1)
    mask = np.max(sliding_window, axis=1) - np.min(sliding_window, axis=1) < number_of_deviations * min_sigma

    return mask


def divide_array(array: np.ndarray, num_sections: int) -> list[slice]:
    section_length = len(array) // num_sections
    remainder = len(array) % num_sections

    slices = []
    start = 0
    for i in range(num_sections):
        if i < remainder:
            end = start + section_length + 1
        else:
            end = start + section_length
        slices.append(slice(start, end))
        start = end

    return slices


if __name__ == "__main__":
    # main_gas()
    main_liq()
    # main_liq_baseline()
