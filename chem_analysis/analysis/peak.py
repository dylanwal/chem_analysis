from __future__ import annotations
import abc
from typing import Protocol, Callable
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
    def __new__(cls, *args, **kwargs):
        if "parent" in kwargs:
            parent = kwargs["parent"]
        else:
            parent = args[0]
        if hasattr(parent, "_" + cls.__name__):
            # intersect class instantiation and redirect it to another variant of peak integration, eg. SEC version
            return super().__new__(parent._PeakIntegration)

        return super().__new__(cls)

    def __init__(self, parent: PeakParent, label=None, id_: int = None):
        self.id_ = id_
        self.parent = parent
        self.label = label

    def to_dict(self) -> OrderedDict:
        attrs = [i for i in self.__dir__() if not i.startswith("_")]

        attrs.sort()
        dict_ = OrderedDict()
        for attr in attrs:
            attr_ = getattr(self, attr)
            if isinstance(attr_, Callable):
                dict_[attr] = attr_()
            else:
                dict_[attr] = attr_

        return dict_


class PeakDiscrete(Peak):
    def __init__(self, parent: PeakParent, index: int, label=None, id_: int = None):
        super().__init__(parent, label, id_)
        self.index = index

    def __str__(self):
        return f"PeakDiscrete(index:{self.index}, x:{self.parent.x[self.index]:.3f}, y:{self.parent.y[self.index]:.3f})"

    def __repr__(self):
        return self.__str__()

    @property
    def value(self) -> float | int:
        return self.parent.y[self.index]


class PeakContinuous(Peak, abc.ABC):
    def __init__(self, parent: PeakParent, properties: PeakProperties = None, label=None, id_: int = None):
        super().__init__(parent, label, id_)
        if properties is None:
            properties = PeakPropertiesData()
        self.properties = properties
        self.properties.set_parent(self)

    def __str__(self):
        return f"PeakContinuous"

    def __repr__(self):
        return self.__str__()

    @property
    @abc.abstractmethod
    def x(self) -> np.ndarray:
        ...

    @property
    @abc.abstractmethod
    def y(self) -> np.ndarray:
        ...


class PeakContinuousData(PeakContinuous):
    def __init__(self,
                 parent: PeakParent,
                 slice_: slice,
                 properties: PeakProperties = None,
                 label=None,
                 id_: int = None
                 ):
        super().__init__(parent, properties, label, id_)
        self.slice_ = slice_

    def __str__(self):
        return f"PeakContinuousData | {self.slice_}"

    @property
    def x(self) -> np.ndarray:
        return self.parent.x[self.slice_]

    @property
    def y(self) -> np.ndarray:
        return self.parent.y[self.slice_]


class PeakContinuousModel(PeakContinuous):
    CUTOFF = 0.001

    def __init__(self,
                 parent: PeakParent,
                 model: Callable,
                 properties: PeakProperties = None,
                 label=None,
                 id_: int = None,
                 slice_: slice = None
                 ):
        super().__init__(parent, properties, label, id_)
        self.model = model
        self.slice_ = slice_

    def __str__(self):
        return f"PeakContinuousModel | {self.model}"

    def _get_slice(self):
        self.slice_ = utils_math.get_slice_by_nearest_y(self.model(self.parent.x), self.CUTOFF)

    @property
    def x(self) -> np.ndarray:
        if self.slice_ is None:
            self._get_slice()

        return self.parent.x[self.slice_]

    @property
    def y(self) -> np.ndarray:
        return self.model(self.x)
