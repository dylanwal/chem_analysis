from typing import Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from src.chem_analysis.analysis.base_obj.signal_ import Signal
from src.chem_analysis.analysis.utils.plot_format import get_plot_color, add_plot_format


class Chromatogram:
    """
    A grouping of Signals.
    """

    __count = 0

    def __init__(self, data: Union[pd.DataFrame, Signal, list[Signal]], name: str = None):
        super().__init__()
        if isinstance(data, pd.DataFrame):
            data = self._load_from_df(data)
        elif isinstance(data, Signal):
            data = {data.y_label if data.y_label is not None else "y_axis": data}
        elif isinstance(data, list):
            if all(isinstance(dat, Signal) for dat in data):
                data = {dat.y_label if dat.y_label is not None else f"y_axis{i}": dat for i, dat in enumerate(data)}
            else:
                raise ValueError("Invalid type in list")

        self.signals = data
        for sig in self.signals:
            setattr(self, sig.name, sig)

        if name is None:
            name = f"Chromat_{Chromatogram.__count}"
            Chromatogram.__count += 1
        self.name = name

    def __repr__(self) -> str:
        text = f"{self.name}: "
        text += "; ".join(self.names)
        return text

    def __iter__(self) -> Signal:
        for sig in self.signals:
            yield sig

    @property
    def names(self):
        return [i.name for i in self.signals]

    @property
    def y_labels(self):
        return [i.y_label for i in self.signals]

    @property
    def x_label(self):
        return self.signals[0].x_label

    @property
    def num_signals(self):
        return len(self.signals)

    def _load_from_df(self, df: pd.DataFrame) -> list[Signal]:
        """ Converts pandas dataframe into traces. """
        signals = []
        for col in df.columns:
            signals.append(Signal(ser=df[col]))

        return signals

    def to_dataframe(self):
        pass

    def baseline(self, **kwargs):
        for sig in self:
            sig.baseline(**kwargs)

    def despike(self, **kwargs):
        for sig in self:
            sig.despike(**kwargs)

    def smooth(self, **kwargs):
        for sig in self:
            sig.smooth(**kwargs)

    def auto_peak_picking(self, **kwargs):
        for sig in self:
            sig.auto_peak_picking(**kwargs)

    def stats(self):
        pass

    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
             op_peaks: bool = True, **kwargs) -> go.Figure:
        if fig is None:
            fig = go.Figure()

        colors = get_plot_color(self.num_signals)

        for sig, color in zip(self, colors):
            kwargs_ = {"color": color}
            if kwargs:
                kwargs_ = {**kwargs_, **kwargs}
            fig = sig.plot(fig, auto_open=False, auto_format=False, op_peaks=op_peaks, **kwargs_)

        if auto_format:
            add_plot_format(fig, self.x_label, "; ".join(self.y_labels))

        if auto_open:
            fig.write_html(f'temp.html', auto_open=True)

        return fig

    def plot_sep_y(self):
        pass


def local_run():
    from scipy.stats import norm
    nx = 1000
    ny = 3
    x = np.linspace(0, nx, nx)
    y = np.empty((ny, nx))
    for i in range(ny):
        rv = norm(loc=nx * np.random.random(1), scale=10)
        y[i, :] = np.linspace(0, nx, nx) + 20 * np.random.random(nx) * np.random.random(1) + 5000 * rv.pdf(x) * \
                  np.random.random(1)

    df = pd.DataFrame(data=y.T, index=x)
    df.columns = ["RI", "UV", "LS"]
    chro = Chromatogram(df)
    chro.baseline(deg=1)
    chro.auto_peak_picking()
    chro.plot()
    print("done")


if __name__ == "__main__":
    local_run()

    # path = r"C:\Users\nicep\Desktop\Reseach_Post\Data\Polyester\DW1_3\SEC\DW1-3-2[DW1-3].csv"
    # df = pd.read_csv(path, header=0, index_col=0)
    # df = df.iloc[:, :10]
    # df.columns = ["LS1", "LS2", "LS3", "LS4", "LS5", "LS6", "LS7", "LS8", "UV", "RI"]
    # df.index.names = ["time (min)"]