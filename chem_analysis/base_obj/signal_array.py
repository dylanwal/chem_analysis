from collections import OrderedDict

import numpy as np

from chem_analysis.utils.sig_fig import apply_sig_figs
from chem_analysis.algorithms.processing.base import Processor
from chem_analysis.base_obj.signal_ import Signal


class SignalArray:
    """
    A grouping of Signals where each Signal occurred at a different time interval.
    """
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None
                 ):
        self.name = name
        self.x_raw = x
        self.y_raw = y
        self.z_raw = z
        self.x_label = x_label if x_label is not None else "x_axis"
        self.y_label = y_label if y_label is not None else "y_axis"
        self.z_label = z_label if z_label is not None else "z_axis"

        self.processor = Processor()
        self._x = None
        self._y = None

        self.signals = []
        # create signals
        for i in self.z.shape[0]:
            sig = Signal(x=self.x, y=z[i, :], var=y[i], x_label=x_label, y_label=z_label)
            sig.id_ = i
            self.signals.append(sig)

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
    def z(self) -> np.ndarray:
        if not self.processor.processed:
            self._z, self._z = self.processor.run(self.x_raw, self.y_raw)

        return self._z

    @property
    def names(self):
        return [i.name for i in self.signals]

    @property
    def number_of_signals(self):
        return len(self.signals)

    def stats(self) -> list[OrderedDict]:
        dicts_ = []
        for sig in self.signals:
            dicts_.append(sig.stats())

        return dicts_

    def print_stats(self, sig_figs: int = 3, output_str: bool = False, **kwargs):
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

        text = tabulate(rows, stats_list[0].keys(), **kwargs)
        if output_str:
            return text

        print(text)
