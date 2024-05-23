import abc
from typing import Protocol, Collection

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
                 ):
        self.atol_x = atol_x

    def evaluate(self, lib_peak: PeakForPickingInterface, signal_peak: PeakContinuous) -> bool:
        if self.atol_x and lib_peak.pos_x - self.atol_x[0] <= signal_peak.max_loc <= lib_peak.pos_x + self.atol_x[1]:
            return True
        return False


class CriteriaAbsoluteRangeAnd(Criteria):
    def __init__(self,
                 atol_x: tuple[float, float] = None,
                 ):
        self.atol_x = atol_x

    def evaluate(self, lib_peak: PeakForPickingInterface, signal_peak: PeakContinuous) -> bool:
        if self.atol_x:
            if lib_peak.pos_x - self.atol_x[0] <= signal_peak.max_loc <= lib_peak.pos_x + self.atol_x[1]:
                return True
        return False


class CriteriaOr(Criteria):
    def __init__(self, criteria: Collection[Criteria]):
        self.criteria = criteria

    def evaluate(self, lib_peak: PeakForPickingInterface, signal_peak: PeakContinuous) -> bool:
        for criteria in self.criteria:
            if criteria.evaluate(lib_peak, signal_peak):
                return True
        return False


class CriteriaAnd(Criteria):
    def __init__(self, criteria: Collection[Criteria]):
        self.criteria = criteria

    def evaluate(self, lib_peak: PeakForPickingInterface, signal_peak: PeakContinuous) -> bool:
        for criteria in self.criteria:
            if not criteria.evaluate(lib_peak, signal_peak):
                return False
        return True
