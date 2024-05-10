
import numpy as np

from chem_analysis.utils.math import map_argmax_to_original
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.processing.weigths.weights import DataWeight
from chem_analysis.analysis.peak_picking.picking_result import ResultPeaks, ResultPeakArray


def find_peak_largest(
        signal: Signal | Signal2D,
        mask: DataWeight = None,
) -> ResultPeaks | ResultPeakArray:
    if isinstance(signal, Signal):
        return find_peak_largest_single(signal, mask)
    elif isinstance(signal, Signal2D):
        results = ResultPeakArray(signal)
        for i in range(len(signal)):
            result = find_peak_largest_single(signal.get_signal(i, processed=True), mask)
            result.signal = signal
            results.add_result(result)
        return results


def find_peak_largest_single(signal: Signal, mask: DataWeight = None) -> ResultPeaks:
    if mask is not None:
        mask = mask.get_mask(signal.x, signal.y)
        y = signal.y[mask]
    else:
        y = signal.y

    indices_of_peaks = np.argmax(y)

    if mask is not None:
        indices_of_peaks = map_argmax_to_original(indices_of_peaks, mask)
    result = ResultPeaks(signal)
    result.indexes = [indices_of_peaks]

    return result
