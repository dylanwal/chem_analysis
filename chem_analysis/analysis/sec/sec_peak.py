from typing import Protocol

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from chem_analysis.analysis.utils.plot_format import add_plot_format, get_multi_y_axis, get_plot_color
from chem_analysis.analysis.base_obj.calibration import Cal
from chem_analysis.analysis.base_obj.peak import Peak
from chem_analysis.analysis.utils import FIGURE_COUNTER


def cal_Mn_D_from_wi(mw_i: np.ndarray, wi: np.ndarray) -> tuple[float, float]:
    """ calculate Mn and D from wi vs MW data (MW goes low to high) """
    data_points = len(mw_i)

    # flip data if giving backwards
    if mw_i[1] > mw_i[-1]:
        mw_i = np.flip(mw_i)
        wi = np.flip(wi)

    wi_d_mi = np.zeros(data_points)
    wi_m_mi = np.zeros(data_points)
    for i in range(data_points):
        if mw_i[i] != 0:
            wi_d_mi[i] = wi[i] / mw_i[i]
        wi_m_mi[i] = wi[i] * mw_i[i]

    mw_n = np.sum(wi) / np.sum(wi_d_mi)
    mw_w = np.sum(wi_m_mi) / np.sum(wi)
    mw_d = mw_w / mw_n
    return mw_n, mw_d


class SECPeakSupports(Protocol):
    name: str
    result: pd.Series
    result_weight: pd.Series
    cal: Cal


class SECPeak(Peak):
    """ SECPeak

    Extends the Peak class to for SEC (size exclusion chromatograph).

    Attributes
    ----------
    mw_i: np.ndarray
        molecular weight of i-mer
    wi: np.ndarray
        weight fraction of i-mer
    xi: np.ndarray
        mole fraction of i-mer
    mw_n: float
        number average molecular weight
    mw_w: float
        weight average molecular weight
    mw_max: float
        molecular weight at peak max
    mw_d: float
        dispersity of molecular weight
    mw_mean: float
        mean molecular weight distribution (same as mw_n)
    mw_std: float
        standard deviation of molecular weight distribution
    mw_skew: float
        skew of molecular weight distribution
        symmetric: -0.5 to 0.5; moderate skew: -1 to -0.5 or 0.5 to 1; high skew: <-1 or >1;
        positive tailing to higher numbers; negative tailing to smaller numbers
    mw_kurtosis: float
        kurtosis of molecular weight distribution (Fisher)  (Warning: highly sensitive to peak bounds)
        negative: flatter peak; positive: sharp peak
    mw_fwhm: float
        full width half max of molecular weight distribution
    mw_asym: float
        asymmetry factor of molecular weight distribution; distance from the center line of the peak to the back
        slope divided by the distance from the center line of the peak to the front slope
        >1 tailing to larger values; <1 tailing to smaller numbers

    """

    def __init__(self, parent: SECPeakSupports, lb_index: int, hb_index: int):
        """

        Parameters
        ----------
        parent: SECPeakSupports
            parent object that the peak is associated with
        lb_index: int
            lower bound index
        hb_index: int
            higher bound index

        """
        self.mw_i = None
        self.wi = None
        self.xi = None

        self.mw_n = None
        self.mw_w = None
        self.mw_max = None
        self.mw_d = None
        self.mw_mean = None
        self.mw_std = None
        self.mw_skew = None
        self.mw_kurtosis = None
        self.mw_fwhm = None
        self.mw_asym = None

        super().__init__(parent, lb_index, hb_index)

    def __repr__(self):
        return f"peak: {self.id_} at {self.max_loc} with Mn: {self.mw_n} (D: {self.mw_d})"

    def calc(self):
        """ Calculates the molecular weight stats for the peak. """
        super().calc()

        self.mw_i = np.flip(self._parent.result_weight.index[self.slice].to_numpy())
        # 'np.flip' so small mw is first avoids having to flip everything in _update below
        self.wi = np.flip(self._parent.result_weight.iloc[self.slice].to_numpy() / \
                  np.trapz(x=self.mw_i, y=self._parent.result_weight.iloc[self.slice].to_numpy()))

        self.mw_n, self.mw_d = cal_Mn_D_from_wi(mw_i=self.mw_i, wi=self.wi)
        self.mw_w = self.mw_n * self.mw_d
        self.mw_max = self._parent.cal(self.max_loc)

        self.xi = self.wi * self.mw_n / self.mw_i

        self.mw_mean = np.trapz(x=self.mw_i, y=self.mw_i*self.xi)
        self.mw_std = np.sqrt(np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 2))
        self.mw_skew = np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 3) / self.mw_std ** 3
        self.mw_kurtosis = (np.trapz(x=self.mw_i, y=self.xi * (self.mw_i - self.mw_mean) ** 4) / self.mw_std ** 4)- 3
        self.mw_fwhm = self.get_fw(x=self.mw_i, y=self.xi, height=0.5)
        self.mw_asym = self.get_asym(x=self.mw_i, y=self.xi, height=0.1)

    def plot_mw(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True, y_label: str = None,
                title: str = None, spread: float = None, **kwargs) -> go.Figure:
        if fig is None:
            fig = go.Figure()

        colors = get_plot_color(2)

        # add main trace
        plot_kwargs_xi = {
            "x": self.mw_i,
            "y": self.xi,
            "mode": 'lines',
            "connectgaps": True,
            "name": "xi",
            "line": dict(color=colors[0]),
            "yaxis": "y1"
        }

        fig.add_trace(go.Scatter(**plot_kwargs_xi))
        plot_kwargs_wi = {
            "x": self.mw_i,
            "y": self.wi,
            "mode": 'lines',
            "connectgaps": True,
            "name": "wi",
            "line": dict(color=colors[1]),
            "yaxis": "y2"
        }

        fig.add_trace(go.Scatter(**plot_kwargs_wi))

        # adding multiple y-axis
        y_axis_labels = ["y1", "y2"]
        if spread is None:
            spread = 0.05 * len(set(y_axis_labels))
        axis_format = get_multi_y_axis(colors, fig, spread)
        fig.update_layout(**axis_format)

        if auto_open:
            global fig_count
            fig.write_html(f'temp{fig_count}.html', auto_open=True)
            fig_count += 1

        return fig

    def stats(self, op_print: bool = True, op_headers: bool = True, window: int = 150, headers: dict = None,
              num_sig_figs: int = 3):
        if headers is None:
            headers = {  # attribute: print
                "id_": "id", "lb_loc": "low bound", "max_loc": "max", "hb_loc": "high bound", "area": "area",
                "mw_n": "mw_n", "mw_w": "mw_w", "mw_d": "mw_d", "mw_max": "mw_max", "mw_std": "mw_std",
                "mw_skew": "mw_skew", "mw_kurtosis": "mw_kurtosis", "mw_asym": "mw_asym"
            }

        return super().stats(op_print, op_headers, window, headers, num_sig_figs)
