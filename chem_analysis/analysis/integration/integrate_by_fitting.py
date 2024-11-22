import logging

import numpy as np

from chem_analysis.analysis.peak_picking.result_picking import ResultPicking, ResultPicking2D
from chem_analysis.analysis.integration.result_integration import ResultIntegration, ResultIntegration2D
from chem_analysis.analysis.peak import PeakContinuous, PeakDiscrete
from chem_analysis.analysis.line_fitting.fitting_normals import integrate_by_fitting_single, n_normals
from chem_analysis.analysis.integration.boundary_detection import rolling_ball_n_points

logger = logging.getLogger(__name__)


def integrate_by_fitting_normal_distribution_single(
        peak: PeakDiscrete,
        num_normals: tuple[int, int] = (1, 3),
        peak_type: type = PeakContinuous,
        id_: int = 0,
) -> PeakContinuous | None:
    lb, ub = rolling_ball_n_points(peak.index, peak.parent.x, peak.parent.y)
    x = peak.parent.x[lb:ub]
    y = peak.parent.y[lb:ub]
    args = integrate_by_fitting_single(x, y, num_normals)

    if args is None:
        # TODO: add checks
        logger.info("peak skipped as fit not successful")
        return None

    for ii in range(int(len(args)/3)):
        y = n_normals(peak.parent.x, *args[3*ii:3*(ii+1)])
        cutoff = 0.001 * np.max(y)
        max_index = np.argmax(y)
        lb_index = np.argmin(np.abs(y[:max_index] - cutoff))
        ub_index = np.argmin(np.abs(y[max_index:] - cutoff)) + max_index

        peak = peak_type(
                parent=peak.parent,
                x=peak.parent.x[lb_index:ub_index],
                y=peak.parent.y[lb_index:ub_index],
                id_=id_
            )
        peak.args = args[3*ii:3*(ii+1)]

    return peak


def integrate_by_fitting_normal_distribution(
    picking_result: ResultPicking,
    num_normals: tuple[int, int] = (1, 3),
):
    result = ResultIntegration(signal=picking_result.signal)

    if len(picking_result.peaks) == 0:
        logger.warning("No peaks to do boundary detection for.")
        return result

    if hasattr(picking_result.signal, "_PeakIntegration"):
        peak_type = picking_result.signal._PeakIntegration
    else:
        peak_type = PeakContinuous

    for i, peak_ in enumerate(picking_result.peaks):
            peak = integrate_by_fitting_normal_distribution_single(peak_, num_normals, peak_type, id_=i)
            if peak is not None:
                result.add_peak(peak)

    return result
