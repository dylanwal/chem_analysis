import abc

import numpy as np

from chem_analysis.analysis.gc_lc_analysis.retention_library import RetentionTimeLibrary


class RetentionTimeMethods(abc.ABC):
    @abc.abstractmethod
    def __call__(self, lib_times: np.ndarray, r_times: np.ndarray) -> np.ndarray:
        ...


class RTMethodATolerance(RetentionTimeMethods):
    def __init__(self, tol: int | float):
        self.tol = tol

    def __call__(self, lib_times: np.ndarray, r_times: np.ndarray) -> np.ndarray:
        mask = -1*np.ones_like(lib_times, dtype=int)
        for i, t in enumerate(r_times):
            index = np.argmin(np.abs(lib_times - t))
            off = np.abs(lib_times[index] - t)
            if off < self.tol:
                mask[i] = index

        return mask


def search_by_retention(
        library: RetentionTimeLibrary,
        retention_times: float | np.ndarray,
        method: RetentionTimeMethods = RTMethodATolerance(0.1),
) -> list:
    mask = method(library.times, retention_times)

    labels = []
    for i, m in enumerate(mask):
        if m:
            labels.append(library.chemicals[i])
        else:
            labels.append(None)
    return labels
