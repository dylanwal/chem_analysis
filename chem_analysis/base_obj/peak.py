from typing import Protocol
from collections import OrderedDict

import numpy as np

import chem_analysis.algorithms.general_math as general_math
from chem_analysis.utils.sig_fig import apply_sig_figs


class PeakSupports(Protocol):
    name: str
    x: np.ndarray
    y: np.ndarray
    y_normalized_by_area: np.ndarray
    y_normalized_by_peak_max: np.ndarray


class PeakID:
    def __init__(self, peak_id: int, signal_id: int):
        self.peak_id = peak_id
        self.signal_id = signal_id

    def __str__(self):
        return f"{self.signal_id}.{self.peak_id}"


class Peak:
    """ Peak

    a peak of a signal

    Attributes
    ----------
    id_: int
        id of peak
    low_bound_index: int
        index of lower bound
    high_bound_index: int
        index of higher bound
    max_index: int
        index of peak maximum
    area: float
        area under the peak
    mean: float
        average value
    std: float
        standard deviation
    skew: float
        skew
        symmetric: -0.5 to 0.5; moderate skew: -1 to -0.5 or 0.5 to 1; high skew: <-1 or >1;
        positive tailing to higher numbers; negative tailing to smaller numbers
    kurtosis: float
        kurtosis (Fisher) (Warning: highly sensitive to peak bounds)
        negative: flatter peak; positive: sharp peak
    full_width_half_max: float
        full width at half maximum
    asymmetry_factor: float
        asymmetry factor; distance from the center line of the peak to the back slope divided by the distance from the
        center line of the peak to the front slope;
        >1 tailing to larger values; <1 tailing to smaller numbers

    """

    __slots__ = ("id_", "parent", "low_bound_index", "high_bound_index", "max_index", "area", "mean", "std",
                 "skew", "kurtosis", "full_width_half_max", "asymmetry_factor")
    _stats = ["id_", "low_bound_location", "max_location", "high_bound_location", "mean", "std", "skew",
              "kurtosis", "full_width_half_max", "asymmetry_factor", "area"]

    def __init__(self, parent: PeakSupports, low_bound_index: int, high_bound_index: int, id_: PeakID = None):
        self.id_ = id_
        self.parent = parent

        self.low_bound_index = low_bound_index
        self.high_bound_index = high_bound_index

        self.max_index = None
        self.area = None
        self.mean = None
        self.std = None
        self.skew = None
        self.kurtosis = None
        self.full_width_half_max = None
        self.asymmetry_factor = None

        self.calculate_parameters()

    def __repr__(self):
        return f"peak: {self.id_} at {self.max_location}"

    @property
    def max_value(self):
        return self.parent.y[self.max_index]

    @property
    def max_location(self):
        return self.parent.x[self.max_index]

    @property
    def low_bound_location(self):
        return self.parent.x[self.low_bound_index]

    @property
    def high_bound_location(self):
        return self.parent.x[self.high_bound_index]

    @property
    def low_bound_value(self):
        return self.parent.y[self.low_bound_index]

    @property
    def high_bound_value(self):
        return self.parent.y[self.high_bound_index]

    @property
    def x(self) -> np.ndarray:
        return self.parent.x[self.low_bound_index:self.high_bound_index]

    @property
    def y(self) -> np.ndarray:
        return self.parent.y[self.low_bound_index:self.high_bound_index]

    def calculate_parameters(self):
        """ Calculates the stats for the peak. """
        x = self.x
        y = self.y
        y_norm = self.parent.y_normalized_by_area[self.low_bound_index:self.high_bound_index]

        self.max_index = np.argmax(y) + self.low_bound_index
        self.mean = general_math.get_mean_of_pdf(x, y_norm=y_norm)
        self.std = general_math.get_standard_deviation_of_pdf(x, y_norm=y_norm, mean=self.mean)
        self.skew = general_math.get_skew_of_pdf(x, y_norm=y_norm, mean=self.mean, standard_deviation=self.std)
        self.kurtosis = general_math.get_kurtosis_of_pdf(x, y_norm=y_norm, mean=self.mean, standard_deviation=self.std)
        self.full_width_half_max = general_math.get_full_width_at_height(x=x, y=y, height=0.5)
        self.asymmetry_factor = general_math.get_asymmetry_factor(x=x, y=y, height=0.1)
        self.area = np.trapz(x=x, y=y)

    def stats(self) -> OrderedDict:
        dict_ = OrderedDict()
        for stat in self._stats:
            dict_[stat] = getattr(self, stat)
        return dict_

    def print_stats(self, sig_figs: int = 3, output_str: bool = False, **kwargs):
        """ Prints stats out for peak. """
        from tabulate import tabulate

        if "tablefmt" not in kwargs:
            kwargs["tablefmt"] = "simple_grid"

        stats = self.stats()
        values = []
        for value in stats.values():
            if isinstance(value, float) or isinstance(value, int):
                value = apply_sig_figs(value, sig_figs)
            values.append(value)

        text = tabulate(values, stats.keys(), **kwargs)
        if output_str:
            return text

        print(text)
