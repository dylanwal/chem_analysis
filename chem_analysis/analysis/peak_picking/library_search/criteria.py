import abc
from typing import Protocol


class PeakInterface(Protocol):
    pos_x: int | float
    pos_y: int | float | None


class Criteria(abc.ABC):
    @abc.abstractmethod
    def evaluate(self, peak: PeakInterface, x: float | int, y: float | int = None) -> bool:
        ...


class CriteriaAbsoluteRangeOr(Criteria):
    def __init__(self,
                 atol_x: tuple[float, float] = None,
                 atol_y: tuple[float, float] = None,
                 ):
        if atol_x is None and atol_y is None:
            raise ValueError(f"'AbsoluteRange': both atol_x and and atol_y can't be 'None'.")
        self.atol_x = atol_x
        self.atol_y = atol_y

    def evaluate(self, peak: PeakInterface, x: float | int, y: float | int = None) -> bool:
        if self.atol_x and peak.pos_x - self.atol_x[0] <= x <= peak.pos_x + self.atol_x[1]:
            return True
        if self.atol_y and y is not None and peak.pos_y - self.atol_y[0] <= y <= peak.pos_y + self.atol_y[1]:
            return True
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

    def evaluate(self, peak: PeakInterface, x: float | int, y: float | int = None) -> bool:
        if self.atol_x:
            if not (peak.pos_x - self.atol_x[0] <= x <= peak.pos_x + self.atol_x[1]):
                return False
        if self.atol_y and y is not None:
            if not (peak.pos_y - self.atol_y[0] <= y <= peak.pos_y + self.atol_y[1]):
                return False
        return True


# TODO: add relative and and or