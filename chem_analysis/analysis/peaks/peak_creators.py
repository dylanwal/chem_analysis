from typing import Sequence, Any

from chem_analysis.utils.math import get_slice
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.peak import PeakData, PeakModel

def peak_from_span(signal: Signal, span: Sequence[int | float], label: Any = None, to_zero: bool = False) -> PeakData:
    """

    Parameters
    ----------
    signal:
    span:
        [left 'x', right 'x']
    to_zero:
        True: integrate down to zero
        False: draw a line between the left and right bound and integrate down to that
    label:
        anything you want to use as a label

    Returns
    -------
    area under the curve
    """
    if isinstance(span, Sequence):
        if len(span) != 2:
            raise ValueError(f"'span' must have length 2. Given: {len(span)}. or be a slice.")
        if span[0] > span[1]:
            raise ValueError(f"'span[0]' must be smaller than 'span[1]'.")
        slice_ = get_slice(signal.x, *span)
    elif isinstance(span, slice):
        slice_ = get_slice(signal.x, span.start, span.stop)
    else:
        raise ValueError(f"'span' must be a slice or Sequence of len(2).")

    return PeakData(
        parent=signal,
        slice_=slice_,
        to_zero=to_zero,
        label=label
    )


from chem_analysis.analysis.line_fitting.peak_models import CurveModel
from chem_analysis.analysis.line_fitting.fitting_main import fitting_simple

def peak_from_model(
        signal: Any,
        span: Sequence[int | float] | slice,
        model:CurveModel,
        fit_span: Sequence[int | float] | slice | None = None,
        label: Any = None,
        to_zero: bool = False
) -> PeakModel:
    if isinstance(span, Sequence):
        if len(span) != 2:
            raise ValueError(f"'span' must have length 2. Given: {len(span)}. or be a slice.")
        slice_ = get_slice(signal.x, *span)
    elif isinstance(span, slice):
        slice_ = get_slice(signal.x, span.start, span.stop)
    else:
        raise ValueError(f"'span' must be a slice or Sequence of len(2).")

    if fit_span is not None:
        if isinstance(fit_span, Sequence):
            if len(fit_span) != 2:
                raise ValueError(f"'fit_span' must have length 2. Given: {len(span)}. or be a slice.")
            slice_fit = get_slice(signal.x, *fit_span)
        elif isinstance(fit_span, slice):
            slice_fit = get_slice(signal.x, fit_span.start, fit_span.stop)
        else:
            raise ValueError(f"'fit_span' must be a slice or Sequence of len(2).")
    else:
        slice_fit = slice_

    parameters = fitting_simple(
        model,
        x=signal.x[slice_fit],
        y=signal.y[slice_fit]
    )
    peak = PeakModel(
        parent=signal,
        model=model,
        slice_=slice_,
        parameters=parameters,
        label=label,
        to_zero=to_zero
    )

    return peak
