from collections import OrderedDict

import numpy as np

import chem_analysis.algorithms.general_math as general_math
from chem_analysis.algorithms.processing.base import Processor
from chem_analysis.utils.sig_fig import apply_sig_figs
from chem_analysis.base_obj.peak import Peak, PeakID


class Signal:
    """ signal

    A signal is any x-y data.

    Attributes
    ----------
    name: str
        Any name the user wants to add.
    x_label: str
        x-axis label
    y_label: str
        y-axis label
    peaks: list[Peak]
        peaks found in signal
    """
    __count = 0

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 name: str = None,
                 raw_data: bool = True,
                 var=None,
                 _parent=None
                 ):
        """

        Parameters
        ----------
        x: np.ndarray
            x data
        y: np.ndarray
            y data
        x_label: str
            x-axis label
        y_label: str
            y-axis label
        name: str
            user defined name

        Notes
        -----
        * Either 'ser' or 'x' and 'y' are required but not both.

        """
        self.x_raw = x
        self.y_raw = y
        self.id_ = Signal.__count
        Signal.__count += 1
        self.name = name if name is not None else f"signal_{self.id_}"
        self.x_label = x_label if x_label is not None else "x_axis"
        self.y_label = y_label if y_label is not None else "y_axis"

        self._peaks = []
        self._peak_counter = 0
        self._parent = _parent
        self.var = var

        self.processor = Processor()
        self._x = None
        self._y = None

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self)})"
        return text

    def __len__(self) -> int:
        return len(self.x)

    @property
    def x(self) -> np.ndarray:
        if not self.processor.processed:
            self._x, self._y = self.processor.run(self.x_raw, self.y_raw)

        return self._x

    @property
    def y(self) -> np.ndarray:
        if not self.processor.processed:
            self._y, self._y = self.processor.run(self.x_raw, self.y_raw)

        return self._y

    @property
    def y_normalized_by_peak_max(self) -> np.ndarray:
        return self.y/np.max(self.y)

    @property
    def y_normalized_by_area(self) -> np.ndarray:
        return general_math.normalize_by_area(self.x, self.y)

    @property
    def peaks(self) -> list[Peak]:
        return self._peaks

    @property
    def tallest_peak(self) -> None | Peak:
        if len(self.peaks) == 0:
            return None

        return max(self.peaks, key=lambda x: x.max_value)

    def add_peak(self, peak: Peak):
        self.peaks.append(peak)
        peak.id_ = PeakID(self.id_, self._peak_counter)
        self._peak_counter += 1

    def delete_peak(self, peak: Peak | int):
        if isinstance(peak, Peak):
            peak = peak.id_
        self.peaks.pop(peak)

        # re-number peaks
        counter = -1
        for peak in self.peaks:
            counter += 1
            peak.id_.signal_id = counter
        self._peak_counter = counter

    def stats(self) -> list[OrderedDict]:
        out = []
        for peak in self.peaks:
            out.append(peak.stats())

        return out

    def print_stats(self, sig_figs: int = 3, output_str: bool = False, **kwargs):
        """ Prints stats out for peak. """
        from tabulate import tabulate

        if "tablefmt" not in kwargs:
            kwargs["tablefmt"] = "simple_grid"

        stats_list = self.stats()
        rows = []
        for stats in stats_list:
            cols = []
            for value in stats.values():
                if isinstance(value, float) or isinstance(value, int):
                    value = apply_sig_figs(value, sig_figs)
                cols.append(value)
            rows.append(cols)

        text = tabulate(rows, stats_list[0].keys(), **kwargs)
        if output_str:
            return text

        print(text)
