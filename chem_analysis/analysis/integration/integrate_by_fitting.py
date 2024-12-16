import logging

import numpy as np

from chem_analysis.analysis.peak_picking.result_picking import ResultPicking, ResultPicking2D
from chem_analysis.analysis.integration.result_integration import ResultIntegration, ResultIntegration2D
from chem_analysis.analysis.peak import PeakContinuousModel, PeakDiscrete
from chem_analysis.analysis.integration.boundary_detection import rolling_ball_n_points
from chem_analysis.analysis.line_fitting.fitting_main import fitting_simple, fitting_adaptive
from chem_analysis.analysis.line_fitting.peak_models import DistributionNormal
from chem_analysis.utils.math import get_slice

logger = logging.getLogger(__name__)


def integrate_by_fitting_normal(
        peak: PeakDiscrete,
        id_: int = 0,
) -> PeakContinuousModel:
    signal = peak.parent
    lb, ub = rolling_ball_n_points(peak.index, signal.x, signal.y)
    slice_ = get_slice(signal.x, lb, ub)
    scale, mean, std = fitting_adaptive(
            model=DistributionNormal,
            x=signal.x[slice_],
            y=signal.y[slice_],
        )

    return PeakContinuousModel(
        parent=signal,
        # slice_=slice_,
        model=DistributionNormal(scale, mean, std),
        label=label,
        id_=id_
    )
