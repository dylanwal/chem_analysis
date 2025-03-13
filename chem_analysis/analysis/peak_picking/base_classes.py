import abc
from typing import Sequence

import numpy as np

from chem_analysis.utils.code_for_subclassing import MixinSubClassList
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.processing.processing_method import Smoothing


class Detector(MixinSubClassList, abc.ABC):
    """ Finds possible peaks. """

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


class Filter(MixinSubClassList, abc.ABC):
    """ Evaluates with a given index pass a Filter. """

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


def find_peaks(
        signal: Signal,
        detectors: Detector | list[Detector],
        filters: Filter | Sequence[Filter] | None = None,
        smoother: Smoothing | None = None
) -> np.ndarray:
    if not isinstance(detectors, list):
        detectors = [detectors]
    if filters is not None and not isinstance(filters, list):
        filters = [filters]

    x, y = np.copy(signal.x), np.copy(signal.y)
    if smoother is not None:
        x, y = smoother.run(x, y)

    indexes = []
    for detector in detectors:
        indexes.append(detector.run_xy(x, y))

    indexes = np.unique(np.concatenate(indexes))

    if filters is not None:
        for filter_ in filters:
            indexes = filter_.run_xy(x, y, indexes)

    return indexes
