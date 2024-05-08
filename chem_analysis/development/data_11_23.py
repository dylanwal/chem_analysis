import glob
import pathlib
from datetime import datetime

import numpy as np
import plotly.graph_objs as go

import chem_analysis.utils.math
from chem_analysis.mass_spec.parsers.agilent_folder import parse_D_folder
from chem_analysis.processing.baseline import bc_polynomial

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
    fig.add_trace(go.Scatter(x=time_normal[1:]/60, y=O2[1:], mode='lines+markers', name="O2"))
    fig.add_trace(go.Scatter(x=time_normal[1:]/60, y=CO2[1:], mode='lines+markers', name="CO2"))
    fig.add_trace(go.Scatter(x=time_normal[1:]/60, y=decane[1:], mode='lines+markers', name="decane"))
    fig.add_trace(go.Scatter(x=[time_o2/60, time_o2/60], y=[0, np.max(O2)], mode='lines', name="O2 start"))
    fig.add_trace(go.Scatter(x=[0, 0], y=[0, np.max(O2)], mode='lines', name="start heat"))
    fig.layout.xaxis.title = "<b>time (min)</b>"
    fig.layout.yaxis.title = "<b>ion count</b>"
    fig.layout.yaxis.range = [0, None]
    fig.show()

    fig = go.Figure(layout={"template": template})
    fig.add_trace(go.Scatter(x=data[15][1]['time'], y=np.sum(data[15][1]['data'], axis=1)))
    fig.show()

#######################################################################################################################


def liquid_single(data):
    ini_dict, ms_dict, fid_dict = data

    ms_data = ms_dict['data']
    ms_time = ms_dict['time']
    ms_sum = np.sum(ms_data, axis=1)

    baseline = bc


def liquid(data):
    return
    O2 = np.empty(len(data))
    CO2 = np.empty_like(O2)
    time = np.empty_like(O2)
    for i, f in enumerate(data):
        time[i], O2[i], CO2[i] = gas_ratio_single(f)
        print(i, time[i], O2[i], CO2[i])


def main_liq():
    folder = r"C:\Users\nicep\Desktop\research_wis\data\11\11_23\gc_ms"
    files = glob.glob(folder + "/" + "DJW-11-23-*min-TMS.D")
    files.sort(key=lambda x: int(x[66:-9]))
    files = [files[2], files[9]]
    data = [parse_D_folder(pathlib.Path(f)) for f in files]
    # liquid(data)

    fig = go.Figure(layout={"template": template})
    fig.add_trace(go.Scatter(x=data[1][1]['time'], y=np.sum(data[1][1]['data'], axis=1)))
    fig.show()


if __name__ == "__main__":
    # main_gas()
    main_liq()
