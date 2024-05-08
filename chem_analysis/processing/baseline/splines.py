import numpy as np
from scipy.interpolate import UnivariateSpline

from chem_analysis.processing.processing_method import Baseline
from chem_analysis.processing.weigths.weights import DataWeight


class Spline(Baseline):
    def __init__(self,
                 degree: int = 3,
                 weights: DataWeight | Iterable[DataWeight] = None,
                 mask: DataWeight | Iterable[DataWeight] = None,
                 non_temporal_processing: bool = False,
                 save_result: bool = False
                 ):
        super().__init__(mask, non_temporal_processing, save_result)
        self.weights = weights
        if not (1 <= degree <= 5):
            raise ValueError('Spline.degree must be 1<=degree<=5')
        self.degree = degree

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        # if len(self.y_sub) == len(y):
            # return self.multiplier * self.y_sub

        return UnivariateSpline(x, y, s=0)

        raise NotImplementedError()  # TODO: x-interpolation
