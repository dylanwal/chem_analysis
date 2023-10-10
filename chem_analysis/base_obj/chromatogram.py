from collections import OrderedDict

import numpy as np

from chem_analysis.utils.sig_fig import apply_sig_figs
from chem_analysis.base_obj.signal_ import Signal


class Chromatogram:
    """
    A grouping of Signals that occur over the same time interval.
    """
    __count = 0

    def __init__(self, data: list[Signal] | np.ndarray, name: str = None):
        if isinstance(data, list):
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
