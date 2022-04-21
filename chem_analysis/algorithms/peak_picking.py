from functools import wraps

from scipy.signal import find_peaks


@wraps(find_peaks)
def scipy_find_peaks(y, *args, **kwargs):
    peaks, _ = find_peaks(y, *args, **kwargs)
    return peaks


