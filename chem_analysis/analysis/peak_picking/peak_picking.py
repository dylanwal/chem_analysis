from functools import wraps
from typing import Protocol

import numpy as np
from scipy.signal import find_peaks

from chem_analysis.utils.general_math import map_argmax_to_original
from chem_analysis.processing.weigths.weights import DataWeight


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

    def __iter__(self):
        return iter(self.indexes)

    def __len__(self):
        return len(self.indexes)

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
            if not (limits[1] <= x <= limits[0]):  # flipped cuz MW goes from big to small
                remove_index.append(i)
        result.indexes = np.delete(result.indexes, remove_index)


@wraps(find_peaks)
def scipy_find_peaks(signal: SignalProtocol, ignore_limits: bool = False, weights: DataWeight = None, **kwargs) \
        -> ResultPeakPicking:
    if weights is not None:
        mask = weights.get_mask(signal.x, signal.y)
        y = signal.y[mask]
    else:
        y = signal.y
    indices_of_peaks, _ = find_peaks(y, **kwargs)
    if weights is not None:
        indices_of_peaks = map_argmax_to_original(indices_of_peaks, mask)
    result = ResultPeakPicking(signal)
    result.indexes = indices_of_peaks

    if not ignore_limits:
        apply_limits(signal, result)

    return result


@wraps(find_peaks)
def max_find_peaks(signal: SignalProtocol, ignore_limits: bool = False, weights: DataWeight = None, **kwargs) \
        -> ResultPeakPicking:
    if weights is not None:
        mask = weights.get_mask(signal.x, signal.y)
        y = signal.y[mask]
    else:
        y = signal.y
    indices_of_peaks = np.argmax(y)
    if weights is not None:
        indices_of_peaks = map_argmax_to_original(indices_of_peaks, mask)
    result = ResultPeakPicking(signal)
    result.indexes = [indices_of_peaks]

    if not ignore_limits:
        apply_limits(signal, result)

    return result