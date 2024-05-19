from typing import Iterable

import numpy as np

from chem_analysis.processing.processing_method import Baseline
from chem_analysis.processing.weigths.weights import DataWeight


class Polynomial(Baseline):
    def __init__(self,
                 degree: int = 1,
                 weights: DataWeight | Iterable[DataWeight] = None,
                 non_temporal_processing: bool = False,
                 save_result: bool = False
                 ):
        super().__init__(non_temporal_processing, save_result)
        self.weights = weights
        self.degree = degree

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        params = np.polyfit(x, y, self.degree, w=self.weights.get_weights(x, y))
        func_baseline = np.poly1d(params)
        return func_baseline(x)
