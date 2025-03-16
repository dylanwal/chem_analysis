from typing import Sequence

import numpy as np

from chem_analysis.analysis.peak import PeakContinuous
from chem_analysis.gc_lc.gc_ms_signal import GCMSSignal
from chem_analysis.mass_spec.ms_signal_2D import MSSignal2D
from chem_analysis.mass_spec.ms_signal import MSSignal


import chem_analysis.utils.math as utils_math


def ms_extract_index(
        signal: GCMSSignal | MSSignal2D,
        index: int | slice | np.ndarray,
) -> MSSignal | list[MSSignal]:
    """

    Parameters
    ----------
    signal:
        signal you want to extract ms from
    index:
        index or slice
        np.ndarray [n,2]

    Returns
    -------

    """
    if isinstance(signal, GCMSSignal):
        signal = signal.ms_raw

    if isinstance(index, int):
        return MSSignal(signal.x, signal.z[index, :])
    elif isinstance(index, slice):
        return MSSignal(signal.x, np.sum(signal.z[index, :], axis=0))
    elif isinstance(index, np.ndarray):
        if index.shape[1] != 2:
            raise ValueError("Shape issue. Format as: [[left, right], [left, right], ...]")
        return [MSSignal(signal.x, np.sum(signal.z[bound[0]:bound[1], :], axis=0)) for bound in index]

    raise ValueError("unsupported index")


def ms_extract_span(
        signal: GCMSSignal | MSSignal2D,
        span: int | float | Sequence[int | float]
) -> MSSignal:
    """

    Parameters
    ----------
    signal:
        signal you want to extract ms from
    span:
        time value or slice for the range you want to extract the ms from

    Returns
    -------

    """
    if span is not None and isinstance(span, Sequence) and len(span) != 2:
        raise ValueError("'ms_extract_span' 'span' argument must be of length 2.")

    if isinstance(span, int) or isinstance(span, float):
        index = np.argmin(abs(signal.y - span))
    elif isinstance(span, Sequence):
        index = utils_math.get_slice(signal.y, span[0], span[1])
    else:
        raise TypeError("'ms_extract_simple' unsupported type for 'span'.")

    return ms_extract_index(signal, index)


def ms_extract_peak(
        peak: PeakContinuous,
) -> MSSignal:
    signal: GCMSSignal = peak.parent
    if not isinstance(signal, GCMSSignal):
        raise ValueError(f"Invalid Peak to extract ms from.\n\tpeak type: {type(peak)}, signal type: {type(signal)}")

    return ms_extract_index(signal, peak.bound_slice)
