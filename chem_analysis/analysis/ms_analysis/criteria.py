import abc
from typing import Protocol

from chem_analysis.analysis.peak import PeakContinuous


class PeakForPickingInterface(Protocol):
    pos_x: int | float
    pos_y: int | float


class Criteria(abc.ABC):
    @abc.abstractmethod
    def evaluate(self, lib_peak: PeakForPickingInterface, signal_peak: PeakContinuous) -> bool:
        ...

#TODO: expand for 2D


class CriteriaAbsoluteRangeOr(Criteria):
    def __init__(self,
                 atol_x: tuple[float, float] = None,
                 atol_y: tuple[float, float] = None,
                 ):
        if atol_x is None and atol_y is None:
            raise ValueError(f"'AbsoluteRange': both atol_x and and atol_y can't be 'None'.")
        self.atol_x = atol_x
        self.atol_y = atol_y

    def evaluate(self, lib_peak: PeakForPickingInterface, signal_peak: PeakContinuous) -> bool:
        if self.atol_x and lib_peak.pos_x - self.atol_x[0] <= signal_peak.max_loc <= lib_peak.pos_x + self.atol_x[1]:
            return True
        # if self.atol_y and y is not None and lib_peak.pos_y - self.atol_y[0] <= y <= lib_peak.pos_y + self.atol_y[1]:
        #     return True
        return False


class CriteriaAbsoluteRangeAnd(Criteria):
    def __init__(self,
                 atol_x: tuple[float, float] = None,
                 atol_y: tuple[float, float] = None,
                 ):
        if atol_x is None and atol_y is None:
            raise ValueError(f"'AbsoluteRange': both atol_x and and atol_y can't be 'None'.")
        self.atol_x = atol_x
        self.atol_y = atol_y

    def evaluate(self, lib_peak: PeakForPickingInterface, signal_peak: PeakContinuous) -> bool:
        if self.atol_x:
            if lib_peak.pos_x - self.atol_x[0] <= signal_peak.max_loc <= lib_peak.pos_x + self.atol_x[1]:
                return True
        # if self.atol_y and y is not None:
        #     if not (lib_peak.pos_y - self.atol_y[0] <= y <= lib_peak.pos_y + self.atol_y[1]):
        #         return False
        return False


# TODO: add relative and and or
