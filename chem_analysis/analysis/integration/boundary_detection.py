import logging

import numpy as np

from chem_analysis.analysis.peak import PeakContinuousData
from chem_analysis.analysis.integration.result_integration import ResultIntegration, ResultIntegration2D

logger = logging.getLogger(__name__)


def rolling_ball_n_points(
        peak_index: int,
        x: np.ndarray,
        y: np.ndarray = None,
        n: int = 2,
        poly_degree: int = 1,
        deriv_degree: int = None,
        max_derivative: float = 0,
        n_points_with_pos_slope: int = 1,
        min_height: float = 0.01,
) -> tuple[int, int]:
    """
    n> 2
    Returns
    -------
    lb_index: int
        index of lower bound
    ub_index: int
        index of upper bound
    """
    min_height = min_height * y[peak_index]
    if deriv_degree is None:
        deriv_degree = poly_degree

    # lower bound
    if peak_index - n <= 0:
        lb_index = peak_index
    else:
        points_with_positive_slope = 0
        for i in range(peak_index - n, 0, -1):
            coefficients = np.polyfit(x[i:i + n], y[i:i + n], poly_degree)
            poly = np.polynomial.Polynomial(list(reversed(coefficients)))
            derivative_value = poly.deriv(deriv_degree)(x[i])
            if derivative_value < max_derivative:  # flip as walking backwards
                points_with_positive_slope += 1
                if points_with_positive_slope >= n_points_with_pos_slope:
                    lb_index = i + points_with_positive_slope
                    break
            else:
                points_with_positive_slope = 0

            if y[i] < min_height:
                lb_index = i
                break
        else:
            lb_index = 0

    # upper bound
    if peak_index + n >= len(x):
        ub_index = peak_index
    else:
        points_with_positive_slope = 0
        for i in range(peak_index + n, len(x)):
            coefficients = np.polyfit(x[i - n:i], y[i - n:i], poly_degree)
            poly = np.polynomial.Polynomial(list(reversed(coefficients)))
            derivative_value = poly.deriv(deriv_degree)(x[i])
            if derivative_value > max_derivative:
                points_with_positive_slope += 1
                if points_with_positive_slope >= n_points_with_pos_slope:
                    ub_index = i - points_with_positive_slope
                    break
            else:
                points_with_positive_slope = 0

            if y[i] < min_height:
                ub_index = i
                break

        else:
            ub_index = len(x) - 1

    return lb_index, ub_index


def rolling_ball(
        peaks,
        n: int = 2,
        poly_degree: int = 1,
        deriv_degree: int = None,
        max_derivative: float = 0,
        n_points_with_pos_slope: int = 1,
        min_height: float = 0.01,
) -> ResultIntegration | ResultIntegration2D:
    """

    Parameters
    ----------
    picking_result
    n
    poly_degree
    deriv_degree
    max_derivative:
        How much it can go up before triggering a bound detection
    n_points_with_pos_slope:
        number of points that can have a slope before triggering
    min_height:
        When to stop if never goes to zero, fraction of max height
    Returns
    -------

    """
    if isinstance(picking_result, ResultPicking):
        return rolling_ball_single(picking_result, n, poly_degree, deriv_degree, max_derivative,
                                   n_points_with_pos_slope, min_height)
    elif isinstance(picking_result, ResultPicking2D):
        result = ResultIntegration2D(signal=picking_result.signal)
        for result_ in picking_result.results:
            result.add_result(
                rolling_ball_single(result_, n, poly_degree, deriv_degree, max_derivative,
                                    n_points_with_pos_slope, min_height)
            )
        return result

    raise TypeError(f"'{rolling_ball.__name__}.picking_result' must be an instance of "
                    f"ResultPicking or ResultPicking2D.")


def rolling_ball_single(
        picking_result,
        n: int,
        poly_degree: int,
        deriv_degree: int,
        max_derivative: float,
        n_points_with_pos_slope: int,
        min_height: float,
):
    result = ResultIntegration(signal=picking_result.signal)

    if len(picking_result.peaks) == 0:
        logger.warning("No peaks to do boundary detection for.")
        return result

    for i, peak in enumerate(picking_result.peaks):
        lb_index, ub_index = rolling_ball_n_points(peak.index, result.signal.x, result.signal.y, n, poly_degree,
                                                   deriv_degree, max_derivative, n_points_with_pos_slope, min_height)

        if not check_end_points(result.signal.y, peak.index, lb_index, ub_index):
            # TODO: improve checks
            continue

        result.add_peak(
            PeakContinuousData(
                parent=picking_result.signal,
                slice_=slice(lb_index, ub_index),
                id_=i
            )
        )
    return result


def check_end_points(y, max_index, lb_index: int, ub_index: int):
    if y[max_index] < y[lb_index] or y[max_index] < y[ub_index]:
        return False
    return True
