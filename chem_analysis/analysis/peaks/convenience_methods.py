from typing import Sequence, Protocol

import numpy as np

from chem_analysis.base_obj.signal_ import Signal


class PeakDetector(Protocol):
    """ Finds peaks. """
    def __call__(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        x:
        y:

        Returns
        -------
        peak index: np.ndarray
            [peak 1, peak 2, ...]
        """


class PeakFilter(Protocol):
    """ Evaluates with a given index pass a PeakFilter. """
    def __call__(self, x: np.ndarray, y: np.ndarray, index: np.ndarray, bounds: np.ndarray | None = None) -> np.ndarray:
        """

        Parameters
        ----------
        x:
        y:
        index: np.ndarray
            index of peak max
        bounds:
            bound of peak
            [[left_index, right_index], [left_index, right_index], ...]

        Returns
        -------
        mask: np.ndarray [bool]
            True: pass filter
            False: don't pass filter

        """


class BoundDetector(Protocol):
    """ Finds edges of peaks. """
    def __call__(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        x:
        y:
        index: np.ndarray
            index of peak max

        Returns
        -------
        bounds: np.ndarray
            [[left_index, right_index], [left_index, right_index], ...]
        """


class Criteria(Protocol):
    def __call__(self, x: np.ndarray, y: np.ndarray, peak: np.ndarray, bound: np.ndarray) -> bool:
        """

        Parameters
        ----------
        x
        y

        Returns
        -------
        bool:
            True: stop
            False: continue, criteria not met

        """


def find_peaks(
        signal: Signal,
        detectors: PeakDetector | list[PeakDetector],
        filters: PeakFilter | Sequence[PeakFilter] | None = None,
) -> np.ndarray:
    return find_peaks_xy(signal.x, signal.y, detectors, filters)


def find_peaks_xy(
        x: np.ndarray,
        y: np.ndarray,
        detectors: PeakDetector | list[PeakDetector],
        filters: PeakFilter | Sequence[PeakFilter] | None = None,
) -> np.ndarray:
    if not isinstance(detectors, list):
        detectors = [detectors]
    if filters is not None and not isinstance(filters, list):
        filters = [filters]

    indices = []
    for detector in detectors:
        indices.append(detector(x, y))

    indices = np.unique(np.concatenate(indices))

    if filters is not None:
        mask = np.ones_like(indices, dtype=bool)
        for filter_ in filters:
            mask &= filter_(x, y, indices)

        indices = indices[mask]

    return indices


def find_peaks_and_bounds(
        signal: Signal,
        peak_det: PeakDetector | list[PeakDetector],
        bound_det: BoundDetector,
        filters: PeakFilter | Sequence[PeakFilter] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    return find_peaks_and_bounds_xy(signal.x, signal.y, peak_det, bound_det, filters)


def find_peaks_and_bounds_xy(
        x: np.ndarray,
        y: np.ndarray,
        peak_det: PeakDetector | list[PeakDetector],
        bound_det: BoundDetector,
        filters: PeakFilter | Sequence[PeakFilter] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    indices = find_peaks_xy(x, y, peak_det, filters)
    bounds = bound_det(x, y, indices)

    if filters is not None:
        mask = np.ones_like(indices, dtype=bool)
        for filter_ in filters:
            mask &= filter_(x, y, indices, bounds)

        indices = indices[mask]

    return indices, bounds


def find_peaks_and_bounds_recursively(
        signal: Signal,
        peak_det: PeakDetector,
        bound_det: BoundDetector,
        stop_criteria: Criteria | Sequence[Criteria] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    return find_peaks_and_bounds_recursively_xy(signal.x, signal.y, peak_det, bound_det, stop_criteria)


def find_peaks_and_bounds_recursively_xy(
        x: np.ndarray,
        y: np.ndarray,
        peak_det: PeakDetector,
        bound_det: BoundDetector,
        stop_criteria: Criteria | Sequence[Criteria],
) -> tuple[np.ndarray, np.ndarray]:
    x_mod, y_mod = np.copy(x), np.copy(y)
    peaks = []
    bounds = []

    for i in range(len(x)):
        peak = get_single_peak(x, y, peak_det(x_mod, y_mod))
        bound = bound_det(x_mod, y_mod, peak)

        for criteria in stop_criteria:
            if criteria(x, y, peak, bound):
                break

        # eliminate peak from y data
        y_mod[bound[0]:bound[1]] = 0

        peaks.append(peak)
        bounds.append(bound)

    return np.array(peaks), np.array(bounds)


def get_single_peak(x: np.ndarray, y: np.ndarray, peaks: np.ndarray) -> np.ndarray:
    """ return one peak"""
    if len(peaks) == 1:
        return peaks
    if len(peaks) == 0:
        return np.argmax(y)

    return peaks[np.argmax(y[peaks])]
