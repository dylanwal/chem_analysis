
import numpy as np


def get_slice(x: np.ndarray, start:  int | float = None, end: int | float = None) -> slice:
    if start is None:
        start_ = 0
    else:
        start_ = np.argmin(np.abs(x-start))
    if end is None:
        end_ = -1
    else:
        end_ = np.argmin(np.abs(x-end))

    return slice(start_, end_)
