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


@wraps(find_peaks)
def scipy_find_peaks(signal: SignalProtocol, *args, **kwargs) -> ResultPeakPicking:
    indices_of_peaks, _ = find_peaks(signal.y, *args, **kwargs)
    result = ResultPeakPicking(signal)
    result.indexes = indices_of_peaks

    return result
