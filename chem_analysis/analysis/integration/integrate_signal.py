

import numpy as np
from scipy.integrate import simpson

from chem_analysis.utils.math import get_slice
from chem_analysis.base_obj.signal_ import Signal


def integrate_trapz(signal: Signal, x_range: tuple[float, float]) -> float | np.ndarray:
    slice_ = get_slice(signal.x, x_range[0], x_range[1])
    return np.trapz(x=signal.x[slice_], y=signal.y[slice_])


def integrate_simpson(signal: Signal, x_range: tuple[float, float]) -> float | np.ndarray:
    slice_ = get_slice(signal.x, x_range[0], x_range[1])
    return simpson(x=signal.x[slice_], y=signal.y[slice_])
