import numpy as np


def poly_baseline(x: np.ndarray, y: np.ndarray, deg: int = 0, mask: np.ndarray = None) \
        -> tuple[np.ndarray, np.ndarray, np.ndarray]:

    if mask is not None:
        x = x[np.where(mask)]
        y = y[np.where(mask)]

    params = np.polyfit(x, y, deg)
    func_baseline = np.poly1d(params)
    y_baseline = func_baseline(x)
    y = y - y_baseline
    return x, y, y_baseline
