import abc
from collections import OrderedDict
import dataclasses
import numpy as np

import chem_analysis.utils.general_math as general_math
from chem_analysis.utils.sig_fig import apply_sig_figs


class Peak(abc.ABC):
    def __init__(self, id_: int = None):
        self.id_ = id_
        self.stats = PeakStats(self)

    @property
    @abc.abstractmethod
    def x(self) -> np.ndarray:
        ...

    @property
    @abc.abstractmethod
    def y(self) -> np.ndarray:
        ...


@dataclasses.dataclass
class PeakParent:
    x: np.ndarray
    y: np.ndarray


class PeakBounded(Peak):
    def __init__(self, parent: PeakParent, bounds: slice, id_: int = None):
        super().__init__(id_)
        self.parent = parent
        self.bounds = bounds

    def __repr__(self):
        return f"peak: {self.id_} at {self.bounds}"

    @property
    def x(self) -> np.ndarray:
        return self.parent.x[self.bounds]

    @property
    def y(self) -> np.ndarray:
        return self.parent.y[self.bounds]

    @property
    def low_bound_value(self) -> float:
        return self.parent.y[self.bounds.start]

    @property
    def high_bound_value(self) -> float:
        return self.parent.y[self.bounds.stop]


class PeakStats:
    """
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
    def __init__(self, parent: Peak):
        self.parent = parent
        self._y_norm = None

    def _get_y_norm(self) -> np.ndarray:
        if self._y_norm is None:
            self._y_norm = self.parent.y/np.trapz(x=self.parent.x, y=self.parent.y)

        return self._y_norm

    @property
    def max_value(self) -> float:
        return np.max(self.parent.y)

    @property
    def max_location(self) -> int:
        return int(np.argmax(self.parent.y))

    @property
    def min_value(self) -> float:
        return np.min(self.parent.y)

    @property
    def min_location(self) -> int:
        return int(np.argmin(self.parent.y))

    @property
    def mean(self) -> float:
        return general_math.get_mean_of_pdf(self.parent.x, y_norm=self._get_y_norm())

    @property
    def std(self):
        return general_math.get_standard_deviation_of_pdf(self.parent.x, y_norm=self._get_y_norm(), mean=self.mean)

    @property
    def skew(self):
        return general_math.get_skew_of_pdf(self.parent.x, y_norm=self._get_y_norm(), mean=self.mean,
                                            standard_deviation=self.std)

    @property
    def kurtosis(self):
        return general_math.get_kurtosis_of_pdf(self.parent.x, y_norm=self._get_y_norm(), mean=self.mean,
                                                standard_deviation=self.std)

    @property
    def full_width_half_max(self):
        return general_math.get_full_width_at_height(x=self.parent.x, y=self.parent.y, height=0.5)

    @property
    def asymmetry_factor(self):
        return general_math.get_asymmetry_factor(x=self.parent.x, y=self.parent.y, height=0.1)

    @property
    def area(self, x: np.ndarray = None) -> float:
        if x is None:
            x = self.parent.x
        return np.trapz(x=x, y=self.parent.y)

    def stats(self) -> OrderedDict:
        dict_ = OrderedDict()
        properties = self.__dict__
        properties.pop("parent")
        properties.pop("_y_norm")
        for stat in self.__dict__:
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
