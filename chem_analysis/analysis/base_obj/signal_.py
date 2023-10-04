from collections import OrderedDict

import numpy as np

import chem_analysis.analysis.algorithms.general_math as general_math
from chem_analysis.analysis.utils.sig_fig import apply_sig_figs
from chem_analysis.analysis.base_obj.peak import Peak


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
    peaks: Peak
        peaks found in signal
    """
    __count = 0
    _peak = Peak

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 name: str = None,
                 raw_data: bool = True,
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
        self.x = x
        self.y = y
        self.id_ = Signal.__count
        Signal.__count += 1
        self.name = name if name is not None else f"signal_{self.id_}"
        self.x_label = x_label if x_label is not None else "x_axis"
        self.y_label = y_label if y_label is not None else "y_axis"

        self.peaks = []
        self._parent = _parent
        self.raw_data = raw_data

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self)})"
        return text

    def __len__(self) -> int:
        return len(self.x)

    @property
    def y_normalized_by_peak_max(self) -> np.ndarray:
        return self.y/np.max(self.y)

    @property
    def y_normalized_by_area(self) -> np.ndarray:
        return general_math.normalize_by_area(self.x, self.y)

    def stats(self) -> list[OrderedDict]:
        out = []
        for peak in self.peaks:
            out.append(peak.stats())
        return out

    def print_stats(self, sig_figs: int = 3, **kwargs):
        """ Prints stats out for peak. """
        from tabulate import tabulate

        if "tablefmt" not in kwargs:
            kwargs["tablefmt"] = "simple_grid"

        stats_list = self.stats()
        rows = []
        for stats in stats_list:
            values = []
            for value in stats.values():
                if isinstance(value, float) or isinstance(value, int):
                    value = apply_sig_figs(value, sig_figs)
                values.append(value)
            rows.append(values)

        print(tabulate(rows, stats_list[0].keys(), **kwargs))
