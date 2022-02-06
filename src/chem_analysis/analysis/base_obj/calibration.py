from typing import Callable

import numpy as np


class Cal:
    def __init__(self, cal: Callable, lb: float = None, ub: float = None, name: str = None):
        self.name = name
        self.cal = cal
        self.lb = lb
        self.ub = ub

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.cal(x)
