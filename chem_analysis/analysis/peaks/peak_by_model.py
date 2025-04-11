from typing import Sequence, Protocol

import numpy as np

from chem_analysis.analysis.peaks.convenience_methods import BoundDetector, Criteria
from chem_analysis.analysis.peaks.bound_detectors import BoundFirstIncreaseZero
from chem_analysis.analysis.line_fitting.fitting_main import fitting_multiple_adaptive, fitting_adaptive, fitting_simple
from chem_analysis.analysis.line_fitting.peak_models import PeakModel, PeakModelBase
import chem_analysis.utils.math as math_utils


def find_peaks_and_bounds_by_model_xy(
        x: np.ndarray,
        y: np.ndarray,
        models: Sequence[PeakModel] | PeakModel,
        stop_criteria: Criteria | Sequence[Criteria],
        bound_det: BoundDetector = None,
        num_trials: int = 3,
        x_range: Sequence[float | int | None, float | int | None] = None
) -> tuple[np.ndarray, np.ndarray]:
    if bound_det is None:
        bound_det = BoundFirstIncreaseZero()

    x_mod, y_mod = np.copy(x), np.copy(y)
    slice_ = math_utils.get_slice(x, *x_range)
    x_mod, y_mod = x_mod[slice_], y_mod[slice_]

    models_ = []
    for i in range(len(x)):
        peak = np.argmax(y)
        bound = bound_det(x_mod, y_mod, peak)

        # fit
        if isinstance(models, Sequence):
            params = fitting_multiple_adaptive(models=models, x=x, y=y, num_trials=num_trials)
        elif isinstance(models, PeakModelBase):
            model = models
            params = fitting_adaptive(model=models, x=x, y=y, num_trials=num_trials)
        else:
            params = fitting_simple(model=models, x=x, y=y, num_trials=num_trials)
            model = models

        for criteria in stop_criteria:
            if criteria(x, y, model, params):
                break

        # eliminate peak from y data
        y_mod[bound[0]:bound[1]] = 0

        models_.append(model)

    return models_
