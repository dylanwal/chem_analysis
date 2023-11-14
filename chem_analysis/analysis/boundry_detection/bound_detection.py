from typing import Protocol
from collections import OrderedDict

import numpy as np

from chem_analysis.analysis.peak_picking.peak_picking import ResultPeakPicking
from chem_analysis.utils.printing_tables import StatsTable


class SignalProtocol(Protocol):
    _peak_type: type
    x: np.ndarray
    y: np.ndarray


class ResultPeakBound:
    def __init__(self, signal: SignalProtocol):
        self.signal = signal
        self.peaks = []

    def __str__(self):
        return f"# of Peaks: {len(self.peaks)}"

    def __repr__(self):
        return self.__str__()

    def get_stats(self) -> list[OrderedDict]:
        dicts_ = []
        for peak in self.peaks:
            dicts_.append(peak.get_stats())

        return dicts_

    def stats_table(self) -> StatsTable:
        return StatsTable.from_list_dicts(self.get_stats())


def rolling_ball_single(
        peak_index: int,
        x: np.ndarray,
        y: np.ndarray = None,
        max_slope: float = 0,
        n_points_with_pos_slope: int = 1,
        min_height: float = 0.01
) -> tuple[int, int]:
    """
    Returns
    -------
    lb_index: int
        index of lower bound
    ub_index: int
        index of upper bound
    """
    min_height = min_height * y[peak_index]
    # lower bound
    if peak_index == 0:
        lb_index = peak_index
    else:
        points_with_positive_slope = 0
        for i in range(peak_index-1, 0, -1):
            slope = (y[i] - y[i + 1]) / (x[i + 1] - x[i])
            if slope > max_slope:
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
    if peak_index == len(x):
        ub_index = peak_index
    else:
        points_with_positive_slope = 0
        for i in range(peak_index + 1, len(x)):
            slope = (y[i] - y[i-1]) / (x[i] - x[i-1])
            if slope > max_slope:
                points_with_positive_slope += 1
                if points_with_positive_slope >= n_points_with_pos_slope:
                    ub_index = i + points_with_positive_slope
                    break
            else:
                points_with_positive_slope = 0

            if y[i] < min_height:
                ub_index = i
                break

        else:
            ub_index = len(x)

    return lb_index, ub_index


def rolling_ball(
        picking_result: ResultPeakPicking,
        n_points_with_pos_slope: int = 1,
        max_slope: float = 0,
        min_height: float = 0.01
) -> ResultPeakBound:
    """
    Parameters
    ----------
    picking_result:

    max_slope: float
        How much it can go up before triggering a bound detection
    n_points_with_pos_slope:
        number of points that can have a slope before triggering
    min_height:
        When to stop if never goes to zero, fraction of max height
    """
    result = ResultPeakBound(picking_result.signal)
    for i, index in enumerate(picking_result.indexes):
        lb_index, ub_index = rolling_ball_single(index, result.signal.x, result.signal.y, max_slope,
                                                 n_points_with_pos_slope, min_height)
        result.peaks.append(
            picking_result.signal._peak_type(
                parent=picking_result.signal,
                bounds=slice(lb_index, ub_index),
                id_=i
            )
        )
    return result
