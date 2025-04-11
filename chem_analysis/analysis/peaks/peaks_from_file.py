import pathlib
import logging

import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.peak import Peak, PeakContinuousData, PeakContinuousModel
from chem_analysis.utils.math import get_slice

logger = logging.getLogger(__name__)


def do_trap(signal: Signal, row_pos: int, params: str) -> PeakContinuousData:
    params = params.split(',')
    if len(params) != 3:
        raise ValueError(f"Row '{row_pos}'| does not have exactly 2 parameters: {len(params)} ({params})")
    try:
        label = params[0]
        lb = float(params[1])
        ub = float(params[2])
    except ValueError:
        raise ValueError(f"Row {row_pos} | parameters must be numbers: {params}")

    slice_ = get_slice(signal.x, lb, ub)
    return PeakContinuousData(
        parent=signal,
        slice_=slice_,
        label=label,
        id_=row_pos
    )


def do_fit_normal(signal: Signal, row_pos: int, params: str) -> PeakContinuousModel:
    params = params.split(',')
    mean, std = None, None
    if not (len(params) == 3 or len(params) == 5):
        raise ValueError(f"Row {row_pos} | 3 or 5 parameters are required: {len(params)} ({params})"
                         f"\n\t Needs at least 'left' and 'right bounds'. ('mean' and 'std' are optional "
                         f"additional params)")
    try:
        label = params[0]
        lb = float(params[1])
        ub = float(params[2])
        if len(params) == 5:
            mean = float(params[3])
            std = float(params[4])
    except ValueError:
        raise ValueError(f"Row {row_pos} | parameters must be numbers: {params}")

    from chem_analysis.analysis.line_fitting.fitting_main import fitting_simple, fitting_adaptive
    from chem_analysis.analysis.line_fitting.peak_models import DistributionNormal
    slice_ = get_slice(signal.x, lb, ub)
    if mean is not None:
        scale, mean, std = fitting_simple(
            model=DistributionNormal,
            x=signal.x[slice_],
            y=signal.y[slice_],
            initial_guess=(np.max(signal.y[slice_]), mean, std)
        )
    else:
        scale, mean, std = fitting_adaptive(
            model=DistributionNormal,
            x=signal.x[slice_],
            y=signal.y[slice_]
        )

    return PeakContinuousModel(
        parent=signal,
        # slice_=slice_,
        model=DistributionNormal(scale, mean, std),
        label=label,
        id_=row_pos
    )


def parse_row(signal: Signal, row_pos: int, row: str) -> Peak:
    row = row.removesuffix('\n')
    try:
        method, params = row.split(',', maxsplit=1)
    except ValueError:
        raise ValueError(f"Row does not have atleast 1 columns (method,label,...): {row}")

    if method == 'trap':
        return do_trap(signal, row_pos, params)  # Validate each row
    elif method == 'fitnormal':
        return do_fit_normal(signal, row_pos, params)

    raise ValueError(f"Row {row_pos} | Invalid 'method': {row}")


def integrate_from_file(
        signal: Signal,
        file_path: str | pathlib.Path,
        raise_errors: bool = True
):
    """

    desecrate peak
    'trap': trapezoidal integration
        provide: left bound [int|float], right bound [int|float]
    'fit_normal': fit normal distribution
        provide: left bound [int|float], right bound [int|float], guess of mean Optional[int|float], guess of std Optional[int|float]
    # 'fitmultinormal': fit multi-normal distribution
    #     provide: left bound [int|float], right bound [int|float], range of normals to consider [int], guess of mean Optional[list[int|float]], guess of std Optional[list[int|float]]
    # 'fitskew': fit skew distribution
    #     provide: left bound [int|float], right bound [int|float], guess of mean Optional[list[int|float]], guess of std Optional[list[int|float]]
    # 'fitdists': fit distributions (normal and skew) distributions
    #         provide: left bound [int|float], right bound [int|float], guess of mean Optional[list[int|float]], guess of std Optional[list[int|float]]
    # 'fitms': fit distributions by deconvolution of ms

    starting rows with '#' can be used it comment-out lines

    Parameters
    ----------
    signal
    file_path
    raise_errors

    Returns
    -------

    """
    result = 1
    errors = []
    with open(file_path, mode='r', encoding='utf-8') as file:
        for i, row in enumerate(file.readlines()):
            if row and not row.startswith("#"):
                try:
                    peak = parse_row(signal, i, row)
                    result.add_peak(peak)
                except Exception as e:
                    logger.error(f"row {i}: {e}")
                    errors.append(e)

    if raise_errors and errors:
        raise errors[0]
    if len(result.peaks) == 0:
        logger.warning("No peaks detection during integration.")

    return result
