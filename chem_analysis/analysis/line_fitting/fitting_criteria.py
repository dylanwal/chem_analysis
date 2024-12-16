from typing import Callable

import numpy as np

Criteria = Callable[[np.ndarray, np.ndarray, np.ndarray, int], int | float]


def BIC(x: np.ndarray, y: np.ndarray, x_model: np.ndarray, num_params: int) -> int | float:
    residuals = y - x_model
    sse = np.sum(residuals ** 2)
    k = num_params
    n = len(x)
    return n * np.log(sse / n) + k * np.log(n)  # BIC formula
