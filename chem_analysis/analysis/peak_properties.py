from __future__ import annotations
from collections import OrderedDict
from typing import Protocol, Callable
from functools import wraps

import numpy as np

import chem_analysis.utils.math as general_math
from chem_analysis.utils.printing_tables import StatsTable


class Parent(Protocol):
    x: np.ndarray
    y: np.ndarray
    
    
class PeakProperties:
    def __init__(self):
        self.parent = None

    def _exclude_from_stats(self):
        return ["_exclude_from_stats", "set_parent", "to_dict", "property_table"]

    def set_parent(self, parent: Parent):
        self.parent = parent

    def to_dict(self) -> OrderedDict:
        attrs = [i for i in self.__dir__() if not i.startswith("_")]

        attrs.sort()
        dict_ = OrderedDict()
        for attr in attrs:
            if attr in self._exclude_from_stats():
                continue
            attr_ = getattr(self, attr)
            if isinstance(attr_, Callable):
                dict_[attr] = attr_()
            else:
                dict_[attr] = attr_

        return dict_

    def property_table(self) -> StatsTable:
        return StatsTable.from_dict(self.to_dict())


class PeakPropertiesData(PeakProperties):
    def __init__(self):
        super().__init__()
        self._y_norm = None

    def _get_y_norm(self) -> np.ndarray:
        if self._y_norm is None:
            self._y_norm = self.parent.y/np.trapz(x=self.parent.x, y=self.parent.y)

        return self._y_norm

    @property
    def min_y(self) -> float:
        return np.min(self.parent.y)

    @property
    def min_index(self) -> int:
        return int(np.argmin(self.parent.y))

    @property
    def min_x(self) -> float:
        return self.parent.x[self.min_index]

    @property
    def max_x(self) -> float:
        return self.parent.x[int(np.argmax(self.parent.y))]

    @property
    def max_index(self) -> int:
        return int(np.argmax(self.parent.y))

    @property
    def max_y(self) -> float:
        return np.max(self.parent.y)

    def area(self, x: np.ndarray = None) -> float:
        if x is None:
            x = self.parent.x
        return np.trapz(x=x, y=self.parent.y)

    def mean(self) -> float:
        return general_math.get_mean_of_pdf(self.parent.x, y_norm=self._get_y_norm())

    def std(self):
        return general_math.get_standard_deviation_of_pdf(self.parent.x, y_norm=self._get_y_norm(), mean=self.mean())

    @wraps(general_math.get_skew_of_pdf)
    def skew(self):
        return general_math.get_skew_of_pdf(self.parent.x, y_norm=self._get_y_norm(), mean=self.mean(),
                                            standard_deviation=self.std())

    @wraps(general_math.get_kurtosis_of_pdf)
    def kurtosis(self):
        return general_math.get_kurtosis_of_pdf(self.parent.x, y_norm=self._get_y_norm(), mean=self.mean(),
                                                standard_deviation=self.std())

    @wraps(general_math.get_full_width_at_height)
    def full_width_half_maximum(self, height: float = 0.5) -> float:
        if not (0 < height < 1):
            raise ValueError('height must be between 0 and 1')
        return general_math.get_full_width_at_height(x=self.parent.x, y=self.parent.y, height=height)

    @wraps(general_math.get_asymmetry_factor)
    def asymmetry_factor(self, height: float = 0.1) -> float:
        if not (0 < height < 1):
            raise ValueError('height must be between 0 and 1')
        return general_math.get_asymmetry_factor(x=self.parent.x, y=self.parent.y, height=height)

    @property
    def bounds(self) -> tuple[float, float]:
        return np.min(self.parent.x), np.max(self.parent.x)

    @property
    def bound_slice(self) -> slice:
        return general_math.get_slice(self.parent.x, self.bounds[0], self.bounds[1])
