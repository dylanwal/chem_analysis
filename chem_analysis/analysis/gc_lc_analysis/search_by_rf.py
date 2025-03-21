import abc

import numpy as np

from chem_analysis.analysis.gc_lc_analysis.retention_library import RetentionTimeLibrary


class RetentionTimeMethods(abc.ABC):
    @abc.abstractmethod
    def __call__(self, lib_times: np.ndarray, r_times: np.ndarray) -> np.ndarray:
        ...


class RTMethodATolerance(RetentionTimeMethods):
    """
    returns the closest compound within tolerance; otherwise it will return None

    """
    def __init__(self, tol: int | float):
        self.tol = tol

    def __call__(self, lib_times: np.ndarray, r_times: np.ndarray) -> np.ndarray:
        indexes = -1*np.ones_like(r_times, dtype=int)
        for i, t in enumerate(r_times):
            index = np.argmin(np.abs(lib_times - t))
            off = np.abs(lib_times[index] - t)
            if off < self.tol:
                indexes[i] = index

        return indexes


class RTMethodNearestN(RetentionTimeMethods):
    """
    returns multiple compounds within tolerance; otherwise it will return None
    """
    def __init__(self, max_matches: int = 2, tol: int | float = 0.1):
        self.max_matches = max_matches
        self.tol = tol

    def __call__(self, lib_times: np.ndarray, r_times: np.ndarray) -> np.ndarray:
        indexes = -1*np.ones((r_times.size, self.max_matches), dtype=int)
        for i, t in enumerate(r_times):
            index = np.argsort(np.abs(lib_times - t))
            off = np.abs(lib_times[index[:self.max_matches]] - t)
            for ii, off_ in enumerate(off):
                if off_ < self.tol:
                    indexes[i, ii] = index[ii]

        return indexes


def search_by_retention(
        library: RetentionTimeLibrary,
        retention_times: float | np.ndarray,
        method: RetentionTimeMethods = RTMethodATolerance(0.1),
) -> list:
    """
    Searches library by retention times.

    Parameters
    ----------
    library:

    retention_times:
        x-values of retention times.
    method:
        method to use for searching.

    Returns
    -------
    list[library compounds] or list[list[library compounds]]
    """
    indexes = method(library.times, retention_times)

    labels = []
    if len(indexes.shape) == 1:
        for i, index in enumerate(indexes):
            if index != -1:
                labels.append(library.chemicals[i])
            else:
                labels.append(None)
    else:
        for i, row in enumerate(indexes):
            sub_labels = []
            for ii, index in enumerate(row):
                if index != -1:
                    sub_labels.append(library.chemicals[index])
            labels.append(sub_labels)

    return labels
