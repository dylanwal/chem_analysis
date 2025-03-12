from typing import Iterable

import numpy as np

from chem_analysis.processing.processing_method import Baseline
from chem_analysis.processing.weigths.weights import Weights


class Polynomial(Baseline):
    def __init__(self,
                 degree: int = 1,
                 weights: Weights | Iterable[Weights] = None,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        super().__init__(temporal_processing, save_result)
        self.weights = weights
        self.degree = degree

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.weights is None:
            params = np.polyfit(x, y, self.degree)
        else:
            params = np.polyfit(x, y, self.degree, w=self.weights.get_weights(x, y))
        func_baseline = np.poly1d(params)
        return func_baseline(x)
