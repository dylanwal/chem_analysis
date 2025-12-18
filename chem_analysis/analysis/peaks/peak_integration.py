from typing import Sequence
import logging

import numpy as np

from chem_analysis.utils.math import get_slice
from chem_analysis.base_obj.signal_ import Signal

logger = logging.getLogger(__name__)

_integral_warning = "Integration includes both positive and negative values which can lead to errors in area calculation."

def base_integrate(x: np.ndarray, y: np.ndarray, to_zero: bool, checks: bool) -> float:
    area = np.trapezoid(x=x, y=y)
    if to_zero:
        if checks and ((y.min() >= 0) or (y.max() <= 0)):
            logger.warning(_integral_warning + f"\n\trange ({x[0]}, {x[-1]})")
    else:
        area_sub = 0.5*(x[-1] - x[0]) *(y[0] + y[-1])
        area -= area_sub
        if checks:
            slope = (y[-1] - y[0]) / (x[-1] - x[0])
            y_baseline = y[0] + slope * (x[-1] - x[0])
            y_corrected = y - y_baseline
            if (y_corrected.min() >= 0) or (y_corrected.max() <= 0):
                logger.warning(_integral_warning + f"\n\trange ({x[0]}, {x[-1]})")

    return area


def integrate_slice_xy(
        x: np.ndarray,
        y: np.ndarray,
        slice_: Sequence[int] | slice,
        to_zero: bool = False,
        checks: bool = True
) -> float | np.ndarray:
    """

    Parameters
    ----------
    x:
    y:
    slice_:
        [left index, right index], slice
    to_zero:
        True: integrate down to zero
        False: draw a line between the left and right bound and integrate down to that
    checks:
        True: checks for positive and negative values in the integral which may cancel each-other out leading to errors.
        False: removes checks; slight performance boost
        provides logger.warning

    Returns
    -------
    area under the curve
    """
    if isinstance(slice_, Sequence) and len(slice_) == 2:
        slice_ = slice(slice_[0], slice_[1])
    x_slice = x[slice_]
    y_slice = y[slice_]
    return base_integrate(x_slice, y_slice, to_zero, checks)


def integrate_slice(
        signal: Signal,
        slice_: np.ndarray | Sequence[int] | slice,
        to_zero: bool = False,
        checks: bool = True
) -> float | np.ndarray:
    """
    Integrate area under the curve given index

    Parameters
    ----------
    signal:
    slice_:
        [left, right] index or slice
    to_zero:
        True: integrate down to zero
        False: draw a line between the left and right bound and integrate down to that
    checks:
        True: checks for positive and negative values in the integral which may cancel each-other out leading to errors.
        False: removes checks; slight performance boost
        provides logger.warning

    Returns
    -------
    area under the curve
    """
    return integrate_slice_xy(signal.x, signal.y, slice_, to_zero, checks)



def integrate(
        signal: Signal,
        span: Sequence[int],
        to_zero: bool = False,
        checks: bool = True
) -> float | np.ndarray:
    """

    Parameters
    ----------
    signal:
    span:
        [left 'x', right 'x']
    to_zero:
        True: integrate down to zero
        False: draw a line between the left and right bound and integrate down to that
    checks:
        True: checks for positive and negative values in the integral which may cancel each-other out leading to errors.
        False: removes checks; slight performance boost
        provides logger.warning

    Returns
    -------
    area under the curve
    """
    return integrate_xy(signal.x, signal.y, span, to_zero, checks)


def integrate_xy(
        x: np.ndarray,
        y: np.ndarray,
        span: Sequence[int],
        to_zero: bool = False,
        checks: bool = True
) -> float | np.ndarray:
    """

    Parameters
    ----------
    x:
    y:
    span:
        [left 'x', right 'x']
    to_zero:
        True: integrate down to zero
        False: draw a line between the left and right bound and integrate down to that
    checks:
        True: checks for positive and negative values in the integral which may cancel each-other out leading to errors.
        False: removes checks; slight performance boost
        provides logger.warning

    Returns
    -------
    area under the curve
    """
    if len(span) != 2:
        raise ValueError(f"'span' must have length 2. Given: {len(span)}")
    slice_ = get_slice(x, *span)
    return integrate_slice_xy(x, y, slice_, to_zero, checks)

