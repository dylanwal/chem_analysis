from typing import Iterable

import numpy as np

from chem_analysis.processing.processing_method import Baseline
from chem_analysis.processing.weigths.weights import DataWeight, DataWeightChain


class BaselineWithMask(Baseline):
    def __init__(self,
                 baseline_method: Baseline,
                 mask: DataWeight | Iterable[DataWeight],
                 non_temporal_processing: bool = False,
                 save_result: bool = False
                 ):
        super().__init__(non_temporal_processing, save_result)
        self.baseline_method = baseline_method
        if mask is not None and isinstance(mask, Iterable):
            mask = DataWeightChain(mask)
        self.mask: DataWeight = mask
        self.mask = mask

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.mask is not None:
            mask = self.mask.get_mask(x, y)
            x_ = x[mask]
            y_ = y[mask]
        else:
            x_ = x
            y_ = y

        baseline = self.baseline_method.get_baseline(x_, y_)
        return np.interp(x, x_, baseline)
