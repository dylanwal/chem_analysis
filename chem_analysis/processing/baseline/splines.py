from typing import Iterable

import numpy as np
from scipy.interpolate import UnivariateSpline

from chem_analysis.processing.processing_method import Baseline
from chem_analysis.processing.weigths.weights import DataWeight


class Spline(Baseline):
    def __init__(self,
                 degree: int = 3,
                 smoothing_factor: float = None,
                 weights: DataWeight | Iterable[DataWeight] = None,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        super().__init__(temporal_processing, save_result)
        self.weights = weights
        if not (1 <= degree <= 5):
            raise ValueError('Spline.degree must be 1<=degree<=5')
        self.degree = degree
        self.smoothing_factor = smoothing_factor

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        spline = UnivariateSpline(x, y, w=self.weights.get_weights(x, y), k=self.degree, s=self.smoothing_factor)
        return spline(x)
