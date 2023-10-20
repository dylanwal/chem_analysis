from typing import Callable
from collections.abc import Sequence

from chem_analysis.base_obj.calibration import Calibration


class SECCalibration(Calibration):
    ...


class ConventionalCalibration(SECCalibration):
    def __init__(self,
                 calibration_function: Callable,
                 *,
                 mw_bounds: Sequence[int | float] = None,
                 time_bounds: Sequence[int | float] = None,
                 name: str = None
                 ):
        super().__init__(calibration_function, y_bounds=mw_bounds, x_bounds=time_bounds, name=name)
