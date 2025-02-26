from __future__ import annotations
import copy

import numpy as np

from chem_analysis.processing.processing_method import ProcessingMethod
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.base_obj.signal_3d import Signal3D


class Processor:
    """
    Processor
    """
    def __init__(self, methods: list[ProcessingMethod] | ProcessingMethod = None):
        if isinstance(methods, ProcessingMethod):
            methods = [methods]
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

    def run(self, sig: Signal | Signal2D | Signal3D) -> Signal | Signal2D | Signal3D:
        new_signal = copy.deepcopy(sig)
        if isinstance(sig, Signal):
            new_signal.x, new_signal.y = self.run_individual(sig.x, sig.y)
        elif isinstance(sig, Signal2D):
            new_signal.x, new_signal.y, new_signal.z = self.run_individual(sig.x, sig.y, sig.z)
        elif isinstance(sig, Signal3D):
            new_signal.x, new_signal.y, new_signal.z, new_signal.w = self.run_individual(sig.x, sig.y, sig.z, sig.w)
        else:
            raise ValueError(f"Unsupported signal type: {type(sig)}")

        return new_signal

    def run_individual(self, x: np.ndarray, y: np.ndarray, z: np.ndarray | None = None, w: np.ndarray | None = None) \
            -> (tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]
                | tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]):
        x = np.copy(x)
        y = np.copy(y)
        if z is not None:
            z = np.copy(z)
        for method in self._methods:
            if z is None:
                x, y = method.run(x, y)
            elif w is None:
                x, y, z = method.run2D(x, y, z)
            else:
                x, y, z, w = method._run3D(x, y, z, w)

        self.processed = True
        if z is None:
            return x, y
        if w is None:
            return x, y, z
        return x, y, z, w

    def get_copy(self) -> Processor:
        copy_ = copy.deepcopy(self)
        copy_.processed = False
        return copy_
