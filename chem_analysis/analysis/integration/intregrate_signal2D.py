import numpy as np
from scipy.integrate import simpson

from chem_analysis.utils.math import get_slice
from chem_analysis.base_obj.signal_2d import Signal2D


def integrate_trapz2D(signal: Signal2D, x_range: tuple[float, float]) -> float | np.ndarray:
    """ along axis 1 """
    slice_ = get_slice(signal.x, x_range[0], x_range[1])
    return np.trapz(x=signal.x[slice_], y=signal.data[:, slice_], axis=1)


def integrate_simpson2D(signal: Signal2D, x_range: tuple[float, float]) -> float | np.ndarray:
    """ along axis 1 """
    slice_ = get_slice(signal.x, x_range[0], x_range[1])
    return simpson(x=signal.x[slice_], y=signal.data[:, slice_], axis=1)
