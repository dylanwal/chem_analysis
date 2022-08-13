from typing import Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from chem_analysis.analysis.base_obj.peak import Peak
from chem_analysis.analysis.utils import logger_analysis, Pipeline, ObjList, FIGURE_COUNTER, up_to_date, array_like
from chem_analysis.analysis.utils.plot_format import get_plot_color, get_similar_color, add_plot_format


class Signal:
    """ signal

    A signal is any x-y data.

    Attributes
    ----------
    name: str
        Any name the user wants to add.
    raw: pd.Series
        raw data
    x_label: str
        x-axis label
    y_label: str
        y-axis label
    pipeline:
        data processing pipeline
    peaks: Peak
        peaks found in signal
    num_peaks: int
        number of peaks
    result: pd.Series
        data process through pipeline
    result_norm: pd.Series
        data process through pipeline normalized (peak max = 1)

    """
    __count = 0
    _peak = Peak

    def __init__(self,
                 name: str = None,
                 ser: pd.Series = None,
                 x: Union[np.ndarray, pd.Series] = None,
                 y: np.ndarray = None,
                 xy: np.ndarray = None,
                 x_label: str = None,
                 y_label: str = None,
                 _parent=None
                 ):
        """

        Parameters
        ----------
        name: str
            user defined name
        ser: pd.Series
            x-y data in the form of a pandas Series
        x: np.ndarray[:]
            x data
        y: np.ndarray[:]
            y data
        xy: np.ndarray[:,2] or [2,:]
            xy data
        x_label: str
            x-axis label
        y_label: str
            y-axis label

        Notes
        -----
        * Either 'ser' or 'x' and 'y' or 'xy' are required but not more than one.

        """
        if name is None:
            name = f"trace_{Signal.__count}"
            Signal.__count += 1

        self.name = name
        self.x_label = x_label
        self.y_label = y_label

        self.raw = self._get_series_from_input(ser, x, y, xy)

        if self.x_label is None:
            self.x_label = "x_axis"
        if self.y_label is None:
            self.y_label = "y_axis"

        self._result = None
        self._result_norm = None
        self._up_to_date = False
        self._parent = _parent
        self.pipeline = Pipeline(up_to_date=self)
        self.peaks = ObjList(self._peak)

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self.raw)})"
        return text

    @property
    @up_to_date
    def result(self) -> pd.Series:
        """ The signal post-processing. """
        return self._result

    @property
    @up_to_date
    def result_norm(self) -> pd.Series:
        """ The signal post-processing normalized with max intensity = 1. """
        self._result_norm = self._result / np.max(self._result)
        return self._result_norm

    @property
    @up_to_date
    def num_peaks(self) -> int:
        return len(self.peaks)

    def _get_series_from_input(self, ser, x, y, xy) -> pd.Series:
        """ Take the data from user and process it into a pd.Series. """
        guard_statements = [
            ser is not None,
            x is not None and y is not None,
            xy is not None]
        if sum(guard_statements) != 1:
            mes = "Too many values provided. Provide either a pandas Series (ser=) or two numpy arrays x and y (x=," \
                  "y=), or a single numpy array (xy=)."
            raise ValueError(f"{mes} (df:{type(ser)}, x:{type(x)}, y:{type(y)}, xy:{type(xy)})")

        # pd.series
        if ser is not None:
            # convert data frame of size 1 to series
            if isinstance(ser, pd.DataFrame) and ser.shape[1] == 1:
                ser = ser[ser.columns[0]]

            if isinstance(ser, pd.Series):
                if self.x_label is None:
                    self.x_label = ser.index.name
                if self.y_label is None:
                    self.y_label = ser.name
                return ser
            else:
                raise TypeError(f"'ser' needs to be a 'pd.Series'. Received '{type(ser)}'.")

        # x,y option
        if x is not None and y is not None:
            if isinstance(x, array_like) and isinstance(y, array_like):
                raw = pd.Series(data=y, index=x, name=self.y_label)
                raw.index.names = [self.x_label]
                return raw
            else:
                raise TypeError(f"'x' and 'y' needs to be a 'np.ndarray'. Received x: '{type(x)}'; y: '{type(y)}'.")

        # xy option
        if xy is not None:
            if isinstance(xy, np.ndarray) and xy.shape[1] == 2:
                raw = pd.Series(data=xy[:, 1], index=xy[:, 0], name=self.y_label)
                raw.index.names = [self.x_label]
                return raw
            else:
                raise TypeError(f"'xy' needs to be a 'np.ndarray'. Received xy: '{type(xy)}'.")

    def _update(self):
        """ performs processing pipeline"""
        x = self.raw.index.to_numpy()
        y = self.raw.to_numpy()
        x, y = self.pipeline.run(x, y)

        self._result = pd.Series(y, x, name=self.y_label)
        self._result.index.names = [self.x_label]
        self._up_to_date = True

    def auto_processing(self):
        import chem_analysis.algorithms.baseline_correction as chem_bc
        self.pipeline.add(chem_bc.adaptive_polynomial_baseline, deg=1)

    @up_to_date
    def peak_picking(self, lb: float, ub: float):
        """ User defines peaks with lower and upper bound. """
        lb_index = int(np.argmin(np.abs(self.raw.index - lb)))
        ub_index = int(np.argmin(np.abs(self.raw.index - ub)))
        self.peaks.add(self._peak(self, lb_index, ub_index))

    @up_to_date
    def auto_peak_picking(self, limit_range: list[float] = None, **kwargs):
        import chem_analysis.algorithms.peak_picking as chem_pp
        import chem_analysis.algorithms.bound_detection as chem_bd
        self.peaks.clear()

        # find peaks
        kwargs_ = {"width": self.raw.index[-1] / 100, "height": 0.03, "prominence": 0.03}
        if kwargs:
            kwargs_ = {**kwargs_, **kwargs}
        if limit_range:
            lb_index = np.argmin(np.abs(self.result.index.to_numpy() - limit_range[0]))
            ub_index = np.argmin(np.abs(self.result.index.to_numpy() - limit_range[1]))
            y = self.result_norm.iloc[lb_index:ub_index].to_numpy()
            y = y / np.max(y)
            peaks_index = chem_pp.scipy_find_peaks(y, **kwargs_) + lb_index
        else:
            peaks_index = chem_pp.scipy_find_peaks(self.result_norm.to_numpy(), **kwargs_)

        # get bounds from peak maximums
        if len(peaks_index) != 0:
            for peak in peaks_index:
                lb, ub = chem_bd.rolling_value(self.result_norm.to_numpy(), peak_index=peak, sensitivity=0.1,
                                               cut_off=0.05)
                self.peaks.add(self._peak(self, lb, ub))
        else:
            logger_analysis.warning(f"No peaks found in signal '{self.name}'.")

        logger_analysis.debug(f"Auto peak picking done on: '{self.name}'. Peaks found: {self.num_peaks}")

    def auto_full(self):
        """
        Does automatic baseline correction and peak detection.
        """
        self.auto_processing()
        self.auto_peak_picking()

    @up_to_date
    def stats(self, op_print: bool = True, op_headers: bool = True, num_sig_figs: int = 3) -> str:
        """ Print out signal/peak stats. """

        text = ""
        for i, peak in enumerate(self.peaks):
            if i == 0:
                if op_headers:
                    text += peak.stats(op_print=False, num_sig_figs=num_sig_figs)
                    continue

            text += peak.stats(op_print=False, op_headers=False, num_sig_figs=num_sig_figs)

        if op_print:
            print(text)
        return text

    @up_to_date
    def stats_df(self) -> pd.DataFrame:
        """ Print out signal/peak stats. """

        series = []
        for peak in self.peaks:
            series.append(peak.stats_series())

        return pd.concat(series, axis=1)

    @up_to_date
    def plot(self, fig: go.Figure = None, auto_open: bool = True, auto_format: bool = True,
             op_peaks: bool = True, y_label: str = None, title: str = None, **kwargs) -> go.Figure:
        """ Plot

        General plotting

        Parameters
        ----------
        fig: go.Figure
            plotly figure; will automatically create if not provided
        auto_open: bool
            create "temp.html" and auto_open in browser
        auto_format: bool
            apply built-in formatting
        op_peaks: bool
            add peak plotting stuff
        y_label: str
            y_axis label (used for multiple y-axis)
        title: str
            title

        Returns
        -------
        fig: go.Figure
            plotly figure

        """
        if fig is None:
            fig = go.Figure()

        if "color" in kwargs:
            color = kwargs.pop("color")
        else:
            color = 'rgb(0,0,0)'

        group = self.name

        # add peaks
        if op_peaks:
            if len(self.peaks) > 0:
                if color == 'rgb(0,0,0)':
                    peak_color = get_plot_color(self.num_peaks)
                else:
                    peak_color = get_similar_color(color, self.num_peaks)

                for peak, color_ in zip(self.peaks, peak_color):
                    peak.plot_add_on(fig, color=color_, group=group, y_label=y_label)

        # add main trace
        plot_kwargs = {
            "x": self.result.index,
            "y": self.result,
            "mode": 'lines',
            "connectgaps": True,
            "name": f"<b>{self.result.name}</b>",
            "legendgroup": group,
            "line": dict(color=color)
        }
        if y_label is not None:
            plot_kwargs["yaxis"] = y_label

        fig.add_trace(go.Scatter(**plot_kwargs))

        if auto_format:
            if title is not None:
                fig.update_layout(title=f"<b>{title}</b>")
            add_plot_format(fig, self.result.index.name, str(self.result.name))

        if auto_open:
            global FIGURE_COUNTER
            fig.write_html(f'temp{FIGURE_COUNTER}.html', auto_open=True)
            FIGURE_COUNTER += 1

        return fig


def local_run():
    from scipy.stats import norm
    n = 1000
    rv = norm(loc=n / 2, scale=10)
    x = np.linspace(0, n, n)
    y = np.linspace(0, n, n) + 20 * np.random.random(n) + 5000 * rv.pdf(x)

    signal = Signal(name="test", x=x, y=y, x_label="x_test", y_label="y_test")
    signal.auto_full()
    signal.stats()
    signal.plot()
    print("done")


if __name__ == '__main__':
    local_run()
