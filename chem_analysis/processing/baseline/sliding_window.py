from typing import Iterable

import numpy as np
from scipy.ndimage import gaussian_filter

from chem_analysis.processing.processing_method import Baseline, Smoothing
from chem_analysis.processing.weigths.weights import DataWeight
from chem_analysis.processing.weigths.sliding_window import sectioned_std


class SectionMinMax(Baseline):
    def __init__(self,
                 window: int = 3,
                 sections: int = 32,
                 number_of_deviations: int | float = 2,
                 smoother: Smoothing = None,
                 mask: DataWeight | Iterable[DataWeight] = None,
                 non_temporal_processing: bool = False,
                 save_result: bool = False
                 ):
        super().__init__(mask, non_temporal_processing, save_result)
        self.window = window
        self.sections = sections
        self.number_of_deviations = number_of_deviations
        self.smoother = smoother

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return baseline_section_std(y, self.window, self.sections, self.number_of_deviations, self.smoother)


def baseline_section_std(y,
                         window: int = 3,
                         sections: int = 32,
                         number_of_deviations: int | float = 2,
                         smoother: Smoothing = None,
                         ):
    mask = sectioned_std(y, window, sections, number_of_deviations, smoother)
    x = np.arange(len(y))
    mask[0], mask[-1] = True, True
    x_mask = x[mask]
    y_mask = y[mask]
    return gaussian_filter(np.interp(x, x_mask, y_mask), 100)

