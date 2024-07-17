from typing import Sequence

import numpy as np


def numpy_to_sparse(array: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    coords = np.flatnonzero(array)
    data = array.ravel()[coords]
    return coords, data, array.shape


def sparse_to_numpy(coords: np.ndarray, data: np.ndarray, shape: Sequence[int]) -> np.ndarray:
    """"""
    x = np.full(shape, 0, data.dtype)
    np.put(x, coords, data)
    return x
