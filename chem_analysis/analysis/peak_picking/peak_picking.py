from functools import wraps
from typing import Protocol

import numpy as np
from scipy.signal import find_peaks


class SignalProtocol(Protocol):
    x: np.ndarray
    y: np.ndarray


class ResultPeakPicking:
    def __init__(self, signal: SignalProtocol):
        self.signal = signal
        self.indexes = None

    def __str__(self):
        return f"# of Peaks: {len(self.indexes)}"

    def __repr__(self):
        return self.__str__()

    def values(self):
        return self.signal.x[self.indexes]

    def stats(self):
        pass


def apply_limits(signal, result: ResultPeakPicking):
    if hasattr(signal, "_limits") and signal._limits() is not None:
        limits = signal._limits()
        remove_index = []
        for i, index in enumerate(result.indexes):
            x = signal.x[index]
            if not (limits[0] <= x <= limits[1]):
                remove_index.append(i)
        result.indexes = np.delete(result.indexes, remove_index)


@wraps(find_peaks)
def scipy_find_peaks(signal: SignalProtocol, ignore_limits: bool = False, **kwargs) -> ResultPeakPicking:
    indices_of_peaks, _ = find_peaks(signal.y, **kwargs)
    result = ResultPeakPicking(signal)
    result.indexes = indices_of_peaks

    if not ignore_limits:
        apply_limits(signal, result)

    return result
