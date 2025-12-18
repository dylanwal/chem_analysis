from typing import Sequence

import numpy as np

from chem_analysis.analysis.peak import PeakData
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
        peak: PeakData,
) -> MSSignal:
    signal: GCMSSignal = peak.parent
    if not isinstance(signal, GCMSSignal):
        raise ValueError(f"Invalid Peak to extract ms from.\n\tpeak type: {type(peak)}, signal type: {type(signal)}")

    return ms_extract_index(signal, peak.bound_slice)



def baseline_corrected_average(ms_data: np.ndarray, rt_slice: slice, mz_slice: slice = slice(None)) -> np.ndarray:
    """
    Extracts and averages a slice of MS data over a retention time range,
    subtracting the baseline defined as values within 3σ of the slice.

    Parameters:
        ms_data (np.ndarray): 2D array representing GC-MS data (n scans x m m/z bins).
        rt_slice (slice): Slice object for selecting retention time range.
        mz_slice (slice): Slice object for selecting m/z range.

    Returns:
        np.ndarray: 1D array of baseline-corrected averaged mass spectrum.
    """
    # Extract the relevant slice
    data_slice = ms_data[rt_slice, mz_slice]

    # Compute baseline as the mean of values within 3σ of the data
    mean_signal = np.mean(data_slice, axis=0)
    std_signal = np.std(data_slice, axis=0)

    # Define baseline as values within 3σ
    mask = (data_slice >= (mean_signal - 3 * std_signal)) & (data_slice <= (mean_signal + 3 * std_signal))

    # Compute baseline as the mean of the masked values
    baseline = np.where(mask, data_slice, np.nan).mean(axis=0)

    # Subtract baseline and compute the final averaged mass spectrum
    corrected_ms = mean_signal - baseline

    return corrected_ms