from typing import Callable

import numpy as np
from scipy.ndimage import gaussian_filter1d

from chem_analysis.processing.processing_method import Baseline
from chem_analysis.processing.smoothing import Smoother
from chem_analysis.processing.weigths.sliding_window_std import sectioned_std

# Input || x: np.ndarray, x_old: np.ndarray, y_old: np.ndarray
# Return || y: np.ndarray
FittingFunction = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]


from scipy.interpolate import CubicSpline


def cubic_splines(x: np.ndarray, x_old: np.ndarray, y_old: np.ndarray, sigma: int = 100) -> np.ndarray:
    y_new = gaussian_filter1d(y_old, sigma)
    cs = CubicSpline(x_old, y_new, extrapolate=True)
    return cs(x)


def interpolation_gaussian_filter(x: np.ndarray, x_old: np.ndarray, y_old: np.ndarray, sigma: int = 100) -> np.ndarray:
    return gaussian_filter1d(np.interp(x, x_old, y_old), sigma)


class SectionMinMax(Baseline):
    def __init__(self,
                 window: int = 3,
                 sections: int = 32,
                 number_of_deviations: int | float = 2,
                 smoother: Smoother = lambda x: gaussian_filter1d(x, 10),
                 fitting_function: FittingFunction = cubic_splines,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        super().__init__(temporal_processing, save_result)
        self.window = window
        self.sections = sections
        self.number_of_deviations = number_of_deviations
        self.smoother = smoother
        self.fitting_function = fitting_function

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return baseline_section_std(y, self.window, self.sections, self.number_of_deviations,
                                    self.smoother, self.fitting_function)


def baseline_section_std(y,
                         window: int = 3,
                         sections: int = 32,
                         number_of_deviations: int | float = 2,
                         smoother: Smoother = lambda x: gaussian_filter1d(x, 10),
                         fitting_function: FittingFunction = cubic_splines,
                         ):
    mask = sectioned_std(y, window, sections, number_of_deviations, smoother)
    x = np.arange(len(y))
    mask[0], mask[-1] = True, True  # include ends
    x_mask = x[mask]
    y_mask = y[mask]
    return fitting_function(x, x_mask, y_mask)



