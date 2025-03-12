import abc
from typing import Sequence

import numpy as np

from chem_analysis.base_obj.signal_ import Signal


class Detector(abc.ABC):
    """ Finds possible peaks. """
    @abc.abstractmethod
    def run(self, x: np.ndarray, y:np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        x
        y

        Returns
        -------
        index of peaks
        """


class Filter(abc.ABC):
    """ Evaluates with a given index pass a Filter. """
    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray, index: np.ndarray) -> np.ndarray:
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
        filters: Filter | Sequence[Filter] | None = None
) -> np.ndarray:
    if not isinstance(detectors, list):
        discovery_methods = list(detectors)
    if filters is not None and not isinstance(filters, list):
        filters = list(filters)

    indexes = []
    for discovery_method in discovery_methods:
        indexes.append(discovery_method.run(signal.x, signal.y))

    indexes = np.array(indexes)

    if filters is not None:
        for filter in filters:
            indexes = filter.run(signal.x, signal.y, indexes)

    return indexes
