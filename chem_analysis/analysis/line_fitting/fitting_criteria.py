from typing import Callable, Protocol, Sequence, Any

import numpy as np

from chem_analysis.analysis.line_fitting.peak_models import CurveModel


class CriteriaFit(Protocol):
    def __call__(self, x: np.ndarray, y: np.ndarray, model: CurveModel, params: Sequence[Any]) -> int | float:
        ...

## goodness of fit



## Choosing between models

def BIC(x: np.ndarray, y: np.ndarray, model: CurveModel, params: Sequence[Any]) -> int | float:
    y_model = model(x, *params)
    residuals = y - y_model
    sse = np.sum(residuals ** 2)
    n = len(x)
    return n * np.log(sse / n) + len(params) * np.log(n)  # BIC formula
