from typing import Callable

from chem_analysis.base_obj.calibration import Calibration


class SECCalibration(Calibration):
    ...


class ConventionalCalibration(SECCalibration):
    def __init__(self,
                 calibration_function: Callable,
                 lower_bound_mw: int | float = None,
                 upper_bound_mw: int | float = None,
                 name: str = None
                 ):
        super().__init__(calibration_function, lower_bound_mw, upper_bound_mw, name)
