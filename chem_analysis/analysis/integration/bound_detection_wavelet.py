import numpy as np
from scipy.signal import find_peaks, cwt, ricker


def detect_peak_boundaries(signal, widths=np.arange(1, 20)):
    # Compute the continuous wavelet transform (CWT) using the Ricker wavelet
    cwt_matrix = cwt(signal, ricker, widths)

    # Use the scale that gives the strongest response
    cwt_sum = np.sum(np.abs(cwt_matrix), axis=0)
    peaks, _ = find_peaks(cwt_sum)

    if len(peaks) == 0:
        return None, None, None

    peak_apex = peaks[np.argmax(signal[peaks])]

    # Find the beginning (left boundary) and end (right boundary)
    left_idx = np.where(signal[:peak_apex] < signal[peak_apex] * 0.1)[0]
    right_idx = np.where(signal[peak_apex:] < signal[peak_apex] * 0.1)[0]

    peak_start = left_idx[-1] if len(left_idx) > 0 else 0
    peak_end = peak_apex + right_idx[0] if len(right_idx) > 0 else len(signal) - 1

    return peak_start, peak_apex, peak_end
