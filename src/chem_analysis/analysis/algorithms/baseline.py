import numpy as np


def poly_baseline(x: np.ndarray, y: np.ndarray, deg: int = 0, mask: np.ndarray = None) \
        -> tuple[np.ndarray, np.ndarray, np.ndarray]:

    if mask is not None:
        x_mask = x[np.where(mask)]
        y_mask = y[np.where(mask)]
    else:
        x_mask = x
        y_mask = y

    params = np.polyfit(x_mask, y_mask, deg)
    func_baseline = np.poly1d(params)
    y_baseline = func_baseline(x)
    y = y - y_baseline
    return x, y, y_baseline


baseline_methods = {"polynomial": poly_baseline}
