import abc
from typing import Sequence, Callable, Any, Optional

import numpy as np

from chem_analysis.utils.code_for_subclassing import MixinSubClassList
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.processing.processing_method import Smoothing

BoundDetector = Callable[[Signal, Optional[Any]], np.ndarray]
BoundDetector_xy = Callable[[np.ndarray, np.ndarray, np.ndarray, Optional[Any]], np.ndarray]


class PeakDetector(MixinSubClassList, abc.ABC):
    """ Finds peaks. """

    def run(self, signal: Signal) -> np.ndarray:
        return self.run_xy(signal.x, signal.y)

    @abc.abstractmethod
    def run_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        x
        y

        Returns
        -------
        index of peaks
        """


class PeakFilter(MixinSubClassList, abc.ABC):
    """ Evaluates with a given index pass a PeakFilter. """

    def run(self, signal: Signal, index: np.ndarray) -> np.ndarray:
        return self.run_xy(signal.x, signal.y, index)

    @abc.abstractmethod
    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        x
        y
        index

        Returns
        -------
        index of peaks
        """


class BoundDetector(MixinSubClassList, abc.ABC):
    """ Finds possible peaks. """

    def run(self, signal: Signal, index: np.ndarray) -> np.ndarray:
        return self.run_xy(signal.x, signal.y, index)

    @abc.abstractmethod
    def run_xy(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        x
        y
        index:
            index of peaks

        Returns
        -------
        index of bounds np.ndarray[n,2]
        """


def find_peaks(
        signal: Signal,
        detectors: PeakDetector | list[PeakDetector],
        filters: PeakFilter | Sequence[PeakFilter] | None = None,
        smoother: Smoothing | None = None
) -> np.ndarray:
    return find_peaks_xy(signal.x, signal.y, detectors, filters, smoother)


def find_peaks_xy(
        x: np.ndarray,
        y: np.ndarray,
        detectors: PeakDetector | list[PeakDetector],
        filters: PeakFilter | Sequence[PeakFilter] | None = None,
        smoother: Smoothing | None = None
) -> np.ndarray:
    if not isinstance(detectors, list):
        detectors = [detectors]
    if filters is not None and not isinstance(filters, list):
        filters = [filters]

    if smoother is not None:
        x, y = np.copy(x), np.copy(y)
        x, y = smoother.run_xy(x, y)

    indexes = []
    for detector in detectors:
        indexes.append(detector.run_xy(x, y))

    indexes = np.unique(np.concatenate(indexes))

    if filters is not None:
        for filter_ in filters:
            indexes = filter_.run_xy(x, y, indexes)

    return indexes


def find_peaks_and_bounds(
        signal: Signal,
        peak: PeakDetector | list[PeakDetector],
        bounds: BoundDetector | BoundDetector_xy,
        filters: PeakFilter | Sequence[PeakFilter] | None = None,
        smoother: Smoothing | None = None
) -> tuple[np.ndarray, np.ndarray]:
    return find_peaks_and_bounds_xy(signal.x, signal.y, peak, bounds, filters, smoother)


def find_peaks_and_bounds_xy(
        x: np.ndarray,
        y: np.ndarray,
        peak: PeakDetector | list[PeakDetector],
        bounds: BoundDetector | BoundDetector_xy,
        filters: PeakFilter | Sequence[PeakFilter] | None = None,
        smoother: Smoothing | None = None
) -> tuple[np.ndarray, np.ndarray]:
    index = find_peaks_xy(x, y, peak, filters, smoother)
    bounds = bounds.run_xy(x, y, index)
    return index, bounds
