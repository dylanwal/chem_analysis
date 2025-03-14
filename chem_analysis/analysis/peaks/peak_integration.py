import numpy as np
from scipy.integrate import simpson

from chem_analysis.utils.math import get_slice
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D


def integrate_trapz(signal: Signal, bounds: np.ndarray) -> float | np.ndarray:
    return integrate_trapz_xy(signal.x, signal.y, bounds)


def integrate_trapz_xy(x: np.ndarray, y: np.ndarray, bounds: np.ndarray) -> float | np.ndarray:
    I = np.empty(bounds.shape[0])
    for i, pair in enumerate(bounds):
        I[i] = np.trapz(x=x[pair[0]:pair[1]], y=y[pair[0]:pair[1]])
    return I


def integrate_trapz_slice(signal: Signal, bounds: np.ndarray) -> float | np.ndarray:
    return integrate_trapz_xy(signal.x, signal.y, bounds)


def integrate_trapz_slice_xy(x: np.ndarray, y: np.ndarray, bounds: np.ndarray) -> float | np.ndarray:
    I = np.empty(bounds.shape[0])
    for i, pair in enumerate(bounds):
        slice_ = get_slice(x, pair[0], pair[1])
        I[i] = np.trapz(x=x[slice_], y=y[slice_])
    return I


def integrate_simpson(signal: Signal, bounds: np.ndarray) -> float | np.ndarray:
    return integrate_simpson_xy(signal.x, signal.y, bounds)


def integrate_simpson_xy(x: np.ndarray, y: np.ndarray, bounds: np.ndarray) -> float | np.ndarray:
    I = np.empty(bounds.shape[0])
    for i, pair in enumerate(bounds):
        I[i] = simpson(x=x[pair[0]:pair[1]], y=y[pair[0]:pair[1]])
    return I


# def integrate_trapz2D(signal: Signal2D, x_range: tuple[float, float]) -> float | np.ndarray:
#     """ along axis 1 """
#     slice_ = get_slice(signal.x, x_range[0], x_range[1])
#     return np.trapz(x=signal.x[slice_], y=signal.data[:, slice_], axis=1)
#
#
# def integrate_simpson2D(signal: Signal2D, x_range: tuple[float, float]) -> float | np.ndarray:
#     """ along axis 1 """
#     slice_ = get_slice(signal.x, x_range[0], x_range[1])
#     return simpson(x=signal.x[slice_], y=signal.data[:, slice_], axis=1)
