from typing import Protocol, Union

import numpy as np
import pandas as pd
import scipy.interpolate as scipy_interpolate
import plotly.graph_objs as go

from src.chem_analysis.analysis.logger import logger_analysis
from src.chem_analysis.analysis.utils.sig_fig import sig_figs


class PeakSupports(Protocol):
    result: pd.Series


def _length_limit(label: str, limit: int) -> str:
    if len(label) > limit:
        return label[0:limit - 5] + "..."

    return label


class Peak:
    def __init__(self, parent: PeakSupports, lb_index: float, hb_index: float, id_: int = None):
        self.id_ = id_
        self._parent = parent

        self.slice_ = slice(lb_index, hb_index)
        self.lb_index = lb_index
        self.hb_index = hb_index
        self.lb_loc = None
        self.hb_loc = None
        self.lb_value = None
        self.hb_value = None

        self.max_index = None
        self.max = None
        self.max_loc = None

        self.area = None
        self.fwhm = None

        self.calc()

    def __repr__(self):
        return f"peak: {self.id_} at {self.max}"

    def calc(self):
        self.lb_value = self._parent.result.iloc[self.lb_index]
        self.hb_value = self._parent.result.iloc[self.hb_index]
        self.lb_loc = self._parent.result.index[self.lb_index]
        self.hb_loc = self._parent.result.index[self.hb_index]

        self.max_index = np.argmax(self._parent.result.iloc[self.slice_]) + self.lb_index
        self.max = self._parent.result.iloc[self.max_index]
        self.max_loc = self._parent.result.index[self.max_index]

        self.area = np.trapz(y=self._parent.result.iloc[self.slice_].to_numpy(),
                             x=self._parent.result.iloc[self.slice_].index.to_numpy())
        self.fwhm = self.get_fwhm(self._parent.result.index.to_numpy(), self._parent.result.to_numpy())

    @staticmethod
    def get_fwhm(x: np.ndarray, y: np.ndarray, k: int = 3) -> Union[None,int]:
        """ Determine full-with-half-maximum of a peaked set of points, x and y. """
        height_half_max = np.max(y) / 2
        s = scipy_interpolate.splrep(x, y - height_half_max, k=k)
        roots = scipy_interpolate.sproot(s)

        if len(roots) > 2:
            msg = "The dataset appears to have multiple peaks, and thus the FWHM can't be determined."
            logger_analysis.warn(msg)
            return None

        elif len(roots) < 2:
            msg = "No proper peaks were found in the data set; likely the dataset is flat (e.g. all zeros)."
            logger_analysis.warn(msg)
            return None

        else:
            return abs(roots[1] - roots[0])

    def _plot(self, fig: go.Figure, group: str = None, **kwargs):
        self._plot_max(fig, group, **kwargs)
        # self._plot_bounds(fig, group, **kwargs)
        self._plot_shade(fig, group, **kwargs)

    def _plot_shade(self, fig: go.Figure, group: str = None, **kwargs):
        kwargs_ = {
            "width": 0
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        if group:
            kkwargs = {"legendgroup": group}
        else:
            kkwargs = {}

        fig.add_trace(go.Scatter(
            x=self._parent.result.index[self.slice_],
            y=self._parent.result.iloc[self.slice_],
            mode="lines",
            fill='tozeroy',
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))

    def _plot_max(self, fig: go.Figure, group: str = None, **kwargs):
        kwargs_ = {
            "size": 1
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        if group:
            kkwargs = {"legendgroup": group}
        else:
            kkwargs = {}

        fig.add_trace(go.Scatter(
            x=[self.max_loc],
            y=[self.max],
            mode="text",
            marker=kwargs_,
            text=[f"peak: {self.id_}"],
            textposition="top center",
            showlegend=False,
            **kkwargs
        ))

    def _plot_bounds(self, fig: go.Figure, group: str = None, height: float = 0.08, **kwargs):
        kwargs_ = {
            "width": 5
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        if group:
            kkwargs = {"legendgroup": group}
        else:
            kkwargs = {}

        bound_height = max(self._parent.result) * height
        # bounds
        fig.add_trace(go.Scatter(
            x=[self.lb_loc, self.lb_loc],
            y=[-bound_height / 2, bound_height / 2],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs

        ))
        fig.add_trace(go.Scatter(
            x=[self.hb_loc, self.hb_loc],
            y=[-bound_height / 2, bound_height / 2],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))
        fig.add_trace(go.Scatter(
            x=[self.lb_loc, self.hb_loc],
            y=[0, 0],
            mode="lines",
            line=kwargs_,
            showlegend=False,
            **kkwargs
        ))

    def stats(self, op_print: bool = True, op_headers: bool = True, window: int = 150):
        text = ""
        headers = {  # attribute: print
            "id_": "id", "lb_loc": "low bound", "max_loc": "max", "hb_loc": "high bound", "area": "area"
        }

        # format
        width = int(window / len(headers))
        row_format = ("{:<" + str(width) + "}") * len(headers)

        # headers
        if op_headers:
            headers_ = [_length_limit(head, width) for head in headers.values()]
            text = row_format.format(*headers_) + "\n"
            text = text + "-" * window + "\n"

        # values
        entries = [_length_limit(str(sig_figs(getattr(self, k), 3)), width) for k in headers]
        text = text + row_format.format(*entries) + "\n"

        if op_print:
            print(text)

        return text
