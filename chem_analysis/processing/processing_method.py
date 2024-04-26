from __future__ import annotations
from typing import Iterable

import abc
import copy

import numpy as np

from chem_analysis.utils.code_for_subclassing import MixinSubClassList


class ProcessingMethod(MixinSubClassList, abc.ABC):

    def __init__(self, non_temporal_processing: bool = False):
        self.non_temporal_processing = non_temporal_processing

    @abc.abstractmethod
    def run(self, x: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...

    def run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.non_temporal_processing:
            return self._run2D(x, y, data)
        return self._run2D_temporal(x, y, data)

    @abc.abstractmethod
    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ...

    def _run2D_temporal(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        for i in range(z.shape[0]):
            _, z[i, :] = self.run(x, z[i, :])
        return x, y, z

    @abc.abstractmethod
    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        ...


class Translation(ProcessingMethod, abc.ABC):
    ...


class Smoothing(ProcessingMethod, abc.ABC):
    ...


class ReSampling(ProcessingMethod, abc.ABC):
    ...


class PhaseCorrection(ProcessingMethod, abc.ABC):
    ...


class FourierTransform(ProcessingMethod, abc.ABC):
    ...


from chem_analysis.processing.weigths.weights import DataWeight, DataWeightChain


class BaselineCorrection(ProcessingMethod, abc.ABC):
    def __init__(self, weights: DataWeight | Iterable[DataWeight] = None, non_temporal_processing: bool = False):
        super().__init__(non_temporal_processing)

        if weights is not None and isinstance(weights, Iterable):
            weights = DataWeightChain(weights)
        self.weights: DataWeight = weights


class Processor:
    """
    Processor
    """
    def __init__(self, methods: list[ProcessingMethod] = None):
        self._methods: list[ProcessingMethod] = [] if methods is None else methods
        self.processed = False

    def __repr__(self):
        return f"Processor: {len(self)} methods"

    def __len__(self):
        return len(self._methods)

    @property
    def methods(self) -> list[ProcessingMethod]:
        return self._methods

    def add(self, *args: ProcessingMethod):
        self._methods += args
        self.processed = False

    def insert(self, index: int, method: ProcessingMethod):
        self._methods.insert(index, method)
        self.processed = False

    def delete(self, method: int | ProcessingMethod):
        if isinstance(method, ProcessingMethod):
            self._methods.remove(method)
        else:
            self._methods.pop(method)
        self.processed = False

    def run(self, x: np.ndarray, y: np.ndarray, z: np.ndarray | None = None) \
            -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]:
        for method in self._methods:
            if z is None:
                x, y = method.run(x, y)
            else:
                x, y, z = method.run2D(x, y, z)

        self.processed = True
        if z is None:
            return x, y
        return x, y, z

    def get_copy(self) -> Processor:
        copy_ = copy.deepcopy(self)
        copy_.processed = False
        return copy_
