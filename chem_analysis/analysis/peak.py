from __future__ import annotations
import abc
from typing import Protocol, Callable, Sequence
from collections import OrderedDict

import numpy as np

import chem_analysis.utils.math as utils_math
from chem_analysis.analysis.peak_properties import PeakProperties, PeakPropertiesData


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
    __slots__ = "id_", "parent", "label"
    def __init__(self, parent: PeakParent, label=None, id_: int = None):
        self.id_ = id_
        self.parent = parent
        self.label = label

    # def to_dict(self) -> OrderedDict:
    #     attrs = [i for i in self.__dir__() if not i.startswith("_")]
    #
    #     attrs.sort()
    #     dict_ = OrderedDict()
    #     for attr in attrs:
    #         attr_ = getattr(self, attr)
    #         if isinstance(attr_, Callable):
    #             dict_[attr] = attr_()
    #         else:
    #             dict_[attr] = attr_
    #
    #     return dict_


class PeakDiscrete(Peak):
    def __init__(self, parent: PeakParent, index: int, label=None, id_: int = None):
        super().__init__(parent, label, id_)
        self.index = index

    def __str__(self):
        return f"{self.__class__.__name__} (index:{self.index}, x:{self.parent.x[self.index]:.3f}, y:{self.parent.y[self.index]:.3f})"

    def __repr__(self):
        return self.__str__()

    @property
    def value(self) -> float | int:
        return self.parent.y[self.index]


class PeakData(Peak):
    __slots__ = "properties", "slice_", "to_zero"
    def __init__(self,
                 parent: PeakParent,
                 slice_: slice,
                 to_zero: bool = False,
                 properties: PeakPropertiesData = None,
                 label=None,
                 id_: int = None
                 ):
        super().__init__(parent, label, id_)
        if properties is None:
                properties = PeakPropertiesData()
        self.properties = properties
        self.properties.set_parent(self)

        self.slice_ = slice_
        self.to_zero = to_zero

    def __str__(self):
        return f"{self.__class__.__name__} | {self.slice_}"

    @property
    def x(self) -> np.ndarray:
        return self.parent.x[self.slice_]

    @property
    def y(self) -> np.ndarray:
        return self.parent.y[self.slice_]

    @property
    def span(self) -> tuple[float | int]:
        return (self.x[self.slice_.start or 0], self.x[self.slice_.stop or -1])


class PeakModel(Peak):
    # CUTOFF = 0.001

    def __init__(self,
                 parent: PeakParent,
                 model: Callable,
                 slice_: slice,
                 parameters: Sequence[int | float],
                 properties: PeakProperties = None,
                 label=None,
                 to_zero: bool = False,
                 id_: int = None
                 ):
        super().__init__(parent, label, id_)
        if properties is None:
            if hasattr(parent, '_peak_properties'):
                properties = getattr(parent, '_peak_properties')()
            else:
                properties = PeakPropertiesData()
        self.properties = properties
        self.properties.set_parent(self)
        self.model = model
        self.slice_ = slice_
        self.parameters = parameters
        self.to_zero = to_zero

    def __str__(self):
        return f"{self.__class__.__name__} | {self.model}"

    @property
    def x(self) -> np.ndarray:
        return self.parent.x[self.slice_]

    @property
    def y(self) -> np.ndarray:
        return self.model(self.x, *self.parameters)

    @property
    def y_parent(self) -> np.ndarray:
        return self.parent.y[self.slice_]

    @property
    def span(self) -> tuple[float | int]:
        return (self.x[self.slice_.start or 0], self.x[self.slice_.stop or -1])
