import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.JSON"
chemistry_lib = ca.mass_spec.GCLibrary.from_JSON(lib_path)

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


def get_time_series(folder_path: str, pattern: str) -> ca.gc_lc.GCMSSignal2D:
    import glob

    files = glob.glob(folder_path + "/" + pattern)
    files.sort(key=lambda x: int(x[66:-9]))
    times = np.array([int(f[66:-9]) for f in files])
    data = [ca.gc_lc.GCParser.from_Agilent_D_folder(f)[0] for f in files]

    return ca.gc_lc.GCMSSignal2D.from_signals(data, times)


def main():
    folder = r"C:\Users\nicep\Desktop\research_wis\data\10\10_12\DJW-10-12-480min-headspace2.D"
    gc_signal, fid_signal = ca.gc_lc.GCParser.from_Agilent_D_folder(folder)

    fig = go.Figure(layout={"template": template})
    fig.add_trace(go.Scatter(x=gc_signal.x, y=gc_signal.y))
    fig.layout.xaxis.title = "<b>time (min)</b>"
    fig.layout.yaxis.title = "<b>ion count</b>"
    fig.show()

    x_slice = ca.utils.math.get_slice(gc_signal.x, 1.78, 1.92)
    fig = go.Figure(layout={"template": template})
    fig.add_trace(go.Bar(x=np.arange(gc_signal.ms_raw.shape[0]), y=np.sum(gc_signal.ms_raw[x_slice, :], axis=0)))
    fig.layout.xaxis.title = "<b>m/z</b>"
    fig.layout.yaxis.title = "<b>ion count</b>"
    fig.layout.xaxis.range = 0, 60
    fig.show()


if __name__ == "__main__":
    main()
