from typing import Sequence

import numpy as np

from chem_analysis.gc_lc.gc_ms_signal import GCMSSignal
from chem_analysis.mass_spec.ms_signal_2D import MSSignal2D
from chem_analysis.mass_spec.ms_signal import MSSignal
from chem_analysis.analysis.integration.result_integration import PeakIntegration


import chem_analysis.utils.math as utils_math


def ms_extract_index(
        signal: GCMSSignal | MSSignal2D,
        index: int | slice,
) -> MSSignal:
    """

    Parameters
    ----------
    signal:
        signal you want to extract ms from
    index:
        index or slice for the ms you want to extract

    Returns
    -------

    """
    if isinstance(signal, GCMSSignal):
        signal = signal.ms_raw

    if isinstance(index, int):
        return MSSignal(signal.x, signal.z[index, :])

    return MSSignal(signal.x, np.sum(signal.z_raw[index, :], axis=0))


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
        peak: PeakIntegration,
) -> MSSignal:
    signal: GCMSSignal = peak.parent
    if not isinstance(signal, GCMSSignal):
        raise ValueError(f"Invalid Peak to extract ms from.\n\tpeak type: {type(peak)}, signal type: {type(signal)}")

    return ms_extract_index(signal, peak.bounds)
