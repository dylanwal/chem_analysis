from typing import Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from chem_analysis.analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.vizilization.plot_format import get_plot_color, add_plot_format, get_multi_y_axis


class Chromatogram:
    """
    A grouping of Signals.
    """

    __count = 0
    _signal = Signal

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
            signals.append(self._signal(ser=df[col], _parent=self))

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

    def auto_peak_baseline(self, **kwargs):
        for sig in self:
            sig.auto_peak_baseline(**kwargs)

    def stats(self, op_print: bool = True) -> str:
        text = ""
        for i, sig in enumerate(self):
            if i == 0:
                text += sig.stats(op_print=False)
                continue

            text += sig.stats(op_print=False, op_headers=False)

        if op_print:
            print(text)
        return text

    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
             op_peaks: bool = True, **kwargs) -> go.Figure:
        if fig is None:
            fig = go.Figure()

        colors = get_plot_color(self.num_signals)

        for sig, color in zip(self, colors):
            kwargs_ = {"color": color}
            if kwargs:
                kwargs_ = {**kwargs_, **kwargs}
            fig = sig.plot_add_on(fig, auto_open=False, auto_format=False, op_peaks=op_peaks, **kwargs_)

        if auto_format:
            add_plot_format(fig, self.x_label, "; ".join(self.y_labels))

        if auto_open:
            fig.write_html(f'temp.html', auto_open=True)

        return fig

    def plot_sep_y(self, fig: go.Figure = None, auto_open: bool = True,
                   op_peaks: bool = True, spread: float = None, **kwargs) -> go.Figure:
        """ Basic plotting """
        if fig is None:
            fig = go.Figure()

        # generate lines
        colors = get_plot_color(self.num_signals)
        y_axis_labels = self._get_y_labels()

        for sig, color, label in zip(self, colors, y_axis_labels):
            kwargs_ = {"color": color}
            if kwargs:
                kwargs_ = {**kwargs_, **kwargs}

            fig = sig.plot_add_on(fig, auto_open=False, auto_format=False, op_peaks=op_peaks, y_label=label, **kwargs_)

        # adding multiple y-axis
        if spread is None:
            spread = 0.05 * len(set(y_axis_labels))
        axis_format = get_multi_y_axis(colors, fig, spread)
        fig.update_layout(**axis_format)

        if auto_open:
            fig.write_html("temp.html", auto_open=True)

        return fig

    def _get_y_labels(self) -> list[str]:
        y_labels = []
        seen = []
        count = 1
        for sig in self:
            if sig.name[:2] in seen:
                index = [i[:2] for i in seen].index(sig.name[:2])
                y_labels.append(y_labels[index])
            else:
                y_labels.append(f"y{count}")
                count += 1

        return y_labels


def local_run():
    from scipy.stats import norm
    nx = 1000
    ny = 3
    x = np.linspace(0, nx, nx)
    y = np.empty((ny, nx))
    for i in range(ny):
        rv = norm(loc=nx * np.random.random(1), scale=10)
        y[i, :] = np.linspace(0, nx, nx) + 20 * np.random.random(nx) * np.random.random(1) + \
                  5000 * rv.pdf(x) * np.random.random(1)

    df = pd.DataFrame(data=y.T, index=x)
    df.columns = ["RI", "UV", "LS"]
    df.index.names = ["time"]
    chro = Chromatogram(df)
    chro.baseline(deg=1)
    chro.auto_peak_picking()
    chro.plot()
    chro.stats()
    print("done")


if __name__ == "__main__":
    local_run()
