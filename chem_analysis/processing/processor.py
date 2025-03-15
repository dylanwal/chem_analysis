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
    def __init__(self, methods: list[ProcessingMethod] | ProcessingMethod = None, *args):
        if isinstance(methods, ProcessingMethod):
            methods = [methods]
        methods = methods or []
        if args:
            methods = methods + list(args)

        self._methods: list[ProcessingMethod] = methods
        self.processed = False

    def __repr__(self):
        return f"Processor: {len(self)} methods"

    def __len__(self):
        return len(self._methods)

    def __copy__(self):
        copy_ = copy.deepcopy(self)
        copy_.processed = False
        return copy_

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
        """ Run processing methods on passed in signal. """
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
        """
        Run processing methods on passed in x, y, (z), (w) np arrays.
        Processing.run_xy(signal) is preferred use over run_individual() method.
        """
        x = np.copy(x)
        y = np.copy(y)
        if z is not None:
            z = np.copy(z)
        if w is not None:
            w = np.copy(w)

        for method in self._methods:
            if z is None:
                x, y = method.run_xy(x, y)
            elif w is None:
                x, y, z = method.run_xyz(x, y, z)
            else:
                x, y, z, w = method._run3D(x, y, z, w)

        self.processed = True

        if z is None:
            return x, y
        elif w is None:
            return x, y, z
        return x, y, z, w
