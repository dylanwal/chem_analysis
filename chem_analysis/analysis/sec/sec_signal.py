from enum import Enum
from typing import Union, Callable

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from chem_analysis.analysis.base_obj.calibration import Calibration
from chem_analysis.analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.sec.sec_peak import SECPeak
from chem_analysis.analysis.utils.plot_format import add_plot_format
from chem_analysis.analysis.utils import FIGURE_COUNTER


class SECType(Enum):
    RI = 0
    UV = 1
    LS = 2
    VISC = 3


class SECSignal(Signal):
    """ SECSignal

    Extends Signal for SEC (size exclusion chronograph)

    Attributes
    ----------
    cal: Calibration
        calibration
    peaks: List[SECPeak]


    """

    _peak = SECPeak

    def __init__(self, cal: Union[Calibration, Callable] = None, **kwrags):
        super().__init__(**kwrags)

        if cal is None and self._parent.cal is not None:
            cal = self._parent.cal
        if not isinstance(cal, Calibration):
            cal = Calibration(cal)
        self.cal = cal

        self._result_weight = None

    @property
    def result_weight(self) -> pd.Series:
        if not self._up_to_date:
            self._update()
        return self._result_weight

    def _update(self):
        super()._update()

        mol_weight = self.cal(self.result.index.to_numpy())
        self._result_weight = pd.Series(self.result.to_numpy(), mol_weight)

    def auto_peak_baseline(self, iterations: int = 3, limit_range: list[float] = None, **kwargs):
        # if limit_range is None and self.cal.lb and self.cal.ub:
        #     limit_range = [self.cal.ub_loc, self.cal.lb_loc]

        super().auto_peak_baseline(iterations=iterations, limit_range=limit_range, **kwargs)

    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True, normalize: bool = False,
             op_peaks: bool = True, op_cal: bool = True, y_label: str = None, title: str = None, **kwargs) -> go.Figure:
        if fig is None:
            fig = go.Figure()

        fig = super().plot(fig, auto_open=False, op_peaks=op_peaks, y_label=y_label, normalize=normalize, **kwargs)

        if op_cal and self.cal is not None:
            self._plot_cal(fig)

        if auto_format:
            if title is not None:
                fig.update_layout(title=title)
            add_plot_format(fig, self.result.index.name, str(self.result.name))
            if op_cal:
                self._add_cal_format(fig)

        if auto_open:
            global FIGURE_COUNTER
            fig.write_html(f'temp{FIGURE_COUNTER}.html', auto_open=True)
            FIGURE_COUNTER += 1

        return fig

    def _plot_cal(self, fig: go.Figure, **kwargs):
        color = "rgb(130,130,130)"
        kwargs_ = {
            "width": 1,
            "color": color
        }
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}

        kkwargs = {}

        mw = self.cal(self.result.index.to_numpy())
        time = self.result.index
        fig.add_trace(go.Scatter(
            x=time,
            y=mw,
            name=f"<b>{self.cal.name if self.cal.name is not None else 'calibration'}</b>",
            mode="lines",
            line=kwargs_,
            yaxis="y2",
            showlegend=True,
            legendgroup="cal",
            **kkwargs
        ))

        kwargs_["dash"] = 'dash'

        # low limit
        if self.cal.lb:
            fig.add_trace(go.Scatter(
                x=[self.cal.lb_loc, self.cal.lb_loc],
                y=[0, np.max(mw)],
                mode="lines",
                line=kwargs_,
                yaxis="y2",
                showlegend=False,
                legendgroup="cal",
                **kkwargs
            ))

        # up limit
        if self.cal.ub:
            fig.add_trace(go.Scatter(
                x=[self.cal.ub_loc, self.cal.ub_loc],
                y=[0, np.max(mw)],
                mode="lines",
                line=kwargs_,
                yaxis="y2",
                showlegend=False,
                legendgroup="cal",
                **kkwargs
            ))

        fig.update_layout(
            yaxis2=dict(
                anchor="x",
                overlaying="y",
                side="right",
                type="log",
                range=[2, 6]
            ),
        )

    def _add_cal_format(self, fig, color: str = 'black'):
        fig.update_layout(
            yaxis2=dict(
                title="<b>molecular weight (g/mol) </b>",
                titlefont=dict(
                    color=color
                ),
                tickfont=dict(
                    color=color
                )
            ),
        )


def local_run():
    from scipy.stats import norm

    def cal_func(time: np.ndarray):
        return 10**(0.0167 * time ** 2 - 0.9225 * time + 14.087)

    cal = Calibration(cal_func, lb=900, ub=319_000)

    nx = 1000
    x = np.linspace(0, 25, nx)
    rv = norm(loc=15, scale=0.6)
    rv2 = norm(loc=18, scale=0.6)
    y = 5 * np.linspace(0, 1, nx) + np.random.random(nx) + 100 * rv.pdf(x) + 20 * rv2.pdf(x)

    sig = SECSignal(name="RI data", x=x, y=y, x_label="time (min)", y_label="intensity", cal=cal)
    sig.auto_peak_baseline(deg=1)
    sig.plot()
    sig.stats()
    print("done")


if __name__ == '__main__':
    local_run()
