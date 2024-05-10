from functools import wraps

import numpy as np
from scipy.signal import find_peaks

from chem_analysis.utils.math import map_argmax_to_original
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.processing.weigths.weights import DataWeight
from chem_analysis.analysis.peak_picking.picking_result import ResultPeakPicking, ResultPeakPickingArray


def apply_mask_signal(signal: Signal, mask: DataWeight = None) -> np.ndarray:
    if mask is not None:
        mask = mask.get_mask(signal.x, signal.y)
        return signal.y[mask]
    return signal.y


def apply_limits(signal, result: ResultPeakPicking):
    if hasattr(signal, "_limits") and signal._limits() is not None:
        limits = signal._limits()
        remove_index = []
        for i, index in enumerate(result.indexes):
            x = signal.x[index]
            if not (limits[1] <= x <= limits[0]):  # flipped cuz MW goes from big to small
                remove_index.append(i)
        result.indexes = np.delete(result.indexes, remove_index)


@wraps(find_peaks)  # TODO: add support for signal array
def find_peaks_scipy(
        signal: Signal | Signal2D,
        mask: DataWeight = None,
        scipy_kwargs: dict = None
) -> ResultPeakPicking | ResultPeakPickingArray:
    if isinstance(signal, Signal):
        return find_peaks_scipy_single(signal, mask, scipy_kwargs)
    elif isinstance(signal, Signal2D):
        results = ResultPeakPickingArray(signal)
        for i in range(len(signal)):
            result = find_peaks_scipy_single(signal.get_signal(i, processed=True), mask, scipy_kwargs)
            # result.signal = signal
            results.add_result(result)
        return results


def find_peaks_scipy_single(signal: Signal, mask: DataWeight = None, scipy_kwargs: dict = None) -> ResultPeakPicking:
    y = apply_mask_signal(signal, mask)

    indices_of_peaks, _ = find_peaks(y, **scipy_kwargs or {})
    if mask is not None:
        indices_of_peaks = map_argmax_to_original(indices_of_peaks, mask)
    result = ResultPeakPicking(signal)
    result.indexes = indices_of_peaks

    return result
