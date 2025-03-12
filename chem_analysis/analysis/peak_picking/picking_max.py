
import numpy as np

from chem_analysis.utils.math import map_argmax_to_original
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.processing.weigths.weights import Weights
from chem_analysis.analysis.peak import PeakDiscrete
from chem_analysis.analysis.peak_picking.result_picking import ResultPicking, ResultPicking2D


def find_peak_largest(
        signal: Signal | Signal2D,
        min_height: float | None = None,
        mask: Weights = None
) -> ResultPicking | ResultPicking2D:
    if isinstance(signal, Signal):
        return find_peak_largest_single(signal, mask, min_height)
    elif isinstance(signal, Signal2D):
        results = ResultPicking2D(signal=signal)
        for sig in signal.signal_iter():
            result = find_peak_largest_single(sig, mask, min_height)
            results.add_result(result)
        return results


def find_peak_largest_single(signal: Signal, mask: Weights = None, min_height: float | None = None) -> ResultPicking:
    result = ResultPicking(signal=signal)
    if mask is not None:
        mask = mask.get_mask(signal.x, signal.y)
        y = signal.y[mask]
    else:
        y = signal.y

    indices_of_peaks = np.argmax(y)

    if mask is not None:
        indices_of_peaks = map_argmax_to_original(indices_of_peaks, mask)

    if min_height is None or signal.y[indices_of_peaks] > min_height:
        result.add_peak(PeakDiscrete(signal, index=indices_of_peaks))
    return result
