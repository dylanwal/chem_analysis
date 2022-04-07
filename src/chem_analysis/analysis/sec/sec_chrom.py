from typing import Union, Callable

import pandas as pd
import plotly.graph_objs as go

from chem_analysis.analysis.base_obj.calibration import Cal
from chem_analysis.analysis.base_obj.chromatogram import Chromatogram
from chem_analysis.analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.sec.sec_signal import SECSignal
from chem_analysis.analysis.utils.plot_format import get_plot_color, add_plot_format
from chem_analysis.analysis.utils import fig_count


class SECChrom(Chromatogram):

    _signal = SECSignal

    def __init__(self, data: Union[pd.DataFrame, Signal, list[Signal]], cal: Union[Cal, Callable] = None):
        if not isinstance(cal, Cal):
            cal = Cal(cal)
        self.cal = cal

        super().__init__(data)

    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
             op_peaks: bool = True, op_cal: bool = True, **kwargs) -> go.Figure:
        if fig is None:
            fig = go.Figure()

        colors = get_plot_color(self.num_signals)

        for i, (sig, color) in enumerate(zip(self, colors)):
            kwargs_ = {"color": color}
            if op_cal and i == 0:
                kwargs_["op_cal"] = True
            else:
                kwargs_["op_cal"] = False
            if kwargs:
                kwargs_ = {**kwargs_, **kwargs}
            fig = sig.plot_add_on(fig, auto_open=False, auto_format=False, op_peaks=op_peaks, **kwargs_)

        if auto_format:
            add_plot_format(fig, self.x_label, "; ".join(self.y_labels))

        if auto_open:
            global fig_count
            fig.write_html(f'temp{fig_count}.html', auto_open=True)
            fig_count += 1

        return fig


def local_run():
    from scipy.stats import norm
    import numpy as np

    def cal_func(time: np.ndarray):
        return 10**(0.0167 * time ** 2 - 0.9225 * time + 14.087)

    cal = Cal(cal_func, lb=900, ub=319_000)

    nx = 1000
    ny = 3
    x = np.linspace(0, 25, nx)
    y = np.empty((ny, nx))
    for i in range(ny):
        rv = norm(loc=15, scale=0.6)
        rv2 = norm(loc=18, scale=0.6)
        y[i, :] = 5 * np.linspace(0, 1, nx) + np.random.random(nx) + 100 * rv.pdf(x) + 20 * rv2.pdf(x)
    df = pd.DataFrame(data=y.T, index=x)
    df.columns = ["RI", "UV", "LS"]
    df.index.names = ["time"]

    chro = SECChrom(data=df, cal=cal)
    chro.auto_peak_baseline(deg=1)
    chro.plot()
    chro.stats()
    print("done")


if __name__ == "__main__":
    local_run()




