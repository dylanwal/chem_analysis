from __future__ import annotations
import abc
import itertools
from collections import OrderedDict
from typing import Protocol, Callable, Iterable
from functools import wraps

import numpy as np

import chem_analysis.utils.math as general_math
from chem_analysis.utils.printing_tables import StatsTable


class PeakParent(Protocol):
    x: np.ndarray
    y: np.ndarray
    name: str


class PeakParent2D(Protocol):
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray


# TODO: add Peak2D give PeakParent2D
class Peak(abc.ABC):
    def __new__(cls, *args, **kwargs):
        if "parent" in kwargs:
            parent = kwargs["parent"]
        else:
            parent = args[0]
        if hasattr(parent, "_" + cls.__name__):
            # intersect class instantiation and redirect it to another variant of peak integration, eg. SEC version
            return super().__new__(parent._PeakIntegration)

        return super().__new__(cls)

    def __init__(self, id_: int = None):
        self.id_ = id_

    def _get_exclude_from_stats(self) -> set[str]:
        return {"stats_table", "stats_dict"}

    def stats_dict(self, exclude: Iterable[str] = None) -> OrderedDict:
        attrs = [i for i in self.__dir__() if not i.startswith("_")]
        if exclude is not None:
            exclude = itertools.chain(self._get_exclude_from_stats(), exclude)
        else:
            exclude = self._get_exclude_from_stats()
        for exclude_ in exclude:
            if exclude_ in attrs:
                attrs.remove(exclude_)

        attrs.sort()
        dict_ = OrderedDict()
        for attr in attrs:
            attr_ = getattr(self, attr)
            if isinstance(attr_, Callable):
                dict_[attr] = attr_()
            else:
                dict_[attr] = attr_

        return dict_

    def stats_table(self, exclude: Iterable[str] = None) -> StatsTable:
        return StatsTable.from_dict(self.stats_dict(exclude))


class PeakDiscrete(Peak):
    def __init__(self, parent: PeakParent, index: int, id_: int = None):
        super().__init__(id_)
        self.parent = parent
        self.index = index

    def _get_exclude_from_stats(self) -> set[str]:
        return super()._get_exclude_from_stats().union({"parent"})

    @property
    def value(self) -> float | int:
        return self.parent.y[self.index]


class PeakContinuous(Peak, abc.ABC):

    def __init__(self, id_: int = None):
        super().__init__(id_)
        self._stats = None
        self._y_norm = None

    @property
    @abc.abstractmethod
    def x(self) -> np.ndarray:
        ...

    @property
    @abc.abstractmethod
    def y(self) -> np.ndarray:
        ...

    def _get_exclude_from_stats(self) -> set[str]:
        return super()._get_exclude_from_stats()

    def _get_y_norm(self) -> np.ndarray:
        if self._y_norm is None:
            self._y_norm = self.y/np.trapz(x=self.x, y=self.y)

        return self._y_norm

    @property
    def min_y(self) -> float:
        return np.min(self.y)

    @property
    def min_index(self) -> int:
        return int(np.argmin(self.y))

    @property
    def min_x(self) -> float:
        return self.x[self.min_index]

    @property
    def max_x(self) -> float:
        return self.x[int(np.argmax(self.y))]

    @property
    def max_index(self) -> int:
        return int(np.argmax(self.y))

    @property
    def max_y(self) -> float:
        return np.max(self.y)

    def area(self, x: np.ndarray = None) -> float:
        if x is None:
            x = self.x
        return np.trapz(x=x, y=self.y)

    def mean(self) -> float:
        return general_math.get_mean_of_pdf(self.x, y_norm=self._get_y_norm())

    def std(self):
        return general_math.get_standard_deviation_of_pdf(self.x, y_norm=self._get_y_norm(), mean=self.mean())

    @wraps(general_math.get_skew_of_pdf)
    def skew(self):
        return general_math.get_skew_of_pdf(self.x, y_norm=self._get_y_norm(), mean=self.mean(),
                                            standard_deviation=self.std())

    @wraps(general_math.get_kurtosis_of_pdf)
    def kurtosis(self):
        return general_math.get_kurtosis_of_pdf(self.x, y_norm=self._get_y_norm(), mean=self.mean(),
                                                standard_deviation=self.std())

    @wraps(general_math.get_full_width_at_height)
    def full_width_half_maximum(self, height: float = 0.5) -> float:
        if not (0 < height < 1):
            raise ValueError('height must be between 0 and 1')
        return general_math.get_full_width_at_height(x=self.x, y=self.y, height=height)

    @wraps(general_math.get_asymmetry_factor)
    def asymmetry_factor(self, height: float = 0.1) -> float:
        if not (0 < height < 1):
            raise ValueError('height must be between 0 and 1')
        return general_math.get_asymmetry_factor(x=self.x, y=self.y, height=height)


class PeakBounded(PeakContinuous):
    def __init__(self, parent: PeakParent, bounds: slice, id_: int = None):
        super().__init__(id_)
        self.parent = parent
        self.bounds = bounds

    def __repr__(self):
        return f"peak: {self.id_} at {self.low_bound_x:.2f}-{self.high_bound_x:.2f}"

    def _get_exclude_from_stats(self) -> set[str]:
        return super()._get_exclude_from_stats().union({"x", "y", "parent"})

    @property
    def x(self) -> np.ndarray:
        return self.parent.x[self.bounds]

    @property
    def y(self) -> np.ndarray:
        return self.parent.y[self.bounds]

    @property
    def low_bound_y(self) -> float:
        return self.parent.y[self.bounds.start]

    @property
    def high_bound_y(self) -> float:
        return self.parent.y[self.bounds.stop]

    @property
    def low_bound_x(self) -> float:
        return self.parent.x[self.bounds.start]

    @property
    def high_bound_x(self) -> float:
        return self.parent.x[self.bounds.stop]

    # @property
    # def area(self) -> float:
    #     return np.trapz(x=self.x, y=self.y)
    #