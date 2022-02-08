from typing import Callable

import numpy as np
import scipy.optimize
import plotly.graph_objs as go


class Cal:
    def __init__(self, cal: Callable, lb: float = None, ub: float = None, name: str = None):
        self.name = name
        self.cal = cal
        self.lb = lb
        self.ub = ub

        self.lb_loc = None
        self.ub_loc = None

        self.calc()

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self.cal(x)

    def calc(self):
        if self.ub:
            def _ub(x: np.ndarray) -> np.ndarray:
                return self(x) - self.ub

            result_ub = scipy.optimize.root_scalar(_ub, x0=0.1, x1=0.2)
            self.ub_loc = result_ub.root

        if self.lb:
            def _lb(x: np.ndarray) -> np.ndarray:
                return self(x) - self.lb

            if self.ub_loc:
                x0 = self.ub_loc
            else:
                x0 = 0.1

            result_lb = scipy.optimize.root_scalar(_lb, x0=x0, x1=x0 + 0.1)
            self.lb_loc = result_lb.root

