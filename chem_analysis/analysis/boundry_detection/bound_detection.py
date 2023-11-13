import abc

import numpy as np

from chem_analysis.analysis import Analysis


class BoundDetection(Analysis, abc.ABC):

    @abc.abstractmethod
    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...


class RollingValue(BoundDetection):
    def __init__(self,
                 peak_index: int = None,
                 peak_value: float | int = None,
                 sensitivity: float = 0.03,
                 cut_off: None | float | int = 0
                 ):
        """
        Parameters
        ----------
        peak_index: int
            index of peak max
        peak_value: int, float
            value of peak max
        sensitivity: float
            How much it can go up before triggering a bound detection, fraction of max height
        cut_off: Union[None, float, int]
            When to stop if never goes to zero, fraction of max height
        """
        self.peak_index = peak_index
        self.peak_value = peak_value
        self.sensitivity = sensitivity
        self.cut_off = cut_off

        if self.peak_index is None:
            raise ValueError("Peak index or peak value are required.")
        if self.peak_value:
            peak_index_ = np.where(x == self.peak_value)
            if len(peak_index_[0]) == 0:
                raise ValueError(f"Peak not found in data. (peak_value: {self.peak_value})")
            elif len(peak_index_[0]) >= 2:
                raise ValueError(f"Multiple peaks have this peak_value. (peak_value: {self.peak_value})")

            peak_index = peak_index_[0]
        else:
            peak_value = x[self.peak_index]

        # set up cutoffs
        if self.cut_off is not None:
            cut_off_value = peak_value * self.cut_off
        else:
            cut_off_value = -10e-50
        if cut_off is not None:
            sensitivity_value = peak_value * self.sensitivity
        else:
            sensitivity_value = 10e-50

    def run(self, x: np.ndarray, y: np.ndarray = None) -> tuple[int, int]:
        """
        Returns
        -------
        lb: int
            index of lower bound
        up: int
            index of upper bound
        """
        # lower bound
        if self.peak_index == 0:
            lb = self.peak_index
        elif self.peak_index < 5:
            lb = np.argmin(x[:self.peak_index])
        else:
            min_ = x[self.peak_index]
            min_index = self.peak_index
            for i, v in enumerate(np.flip(x[:self.peak_index])):
                if v < min_:
                    min_ = v
                    min_index = self.peak_index - i
                    if min_ < self.cut_off_value:
                        break
                if v - x[min_index] > self.sensitivity_value:
                    break

            lb = min_index - 1

        # upper bound
        if self.peak_index == len(x):
            ub = self.peak_index
        elif len(x) - self.peak_index < 5:
            ub = np.argmin(x[self.peak_index:])
        else:
            min_ = x[self.peak_index]
            min_index = self.peak_index
            for i, v in enumerate(x[self.peak_index:]):
                if v < min_:
                    min_ = v
                    min_index = self.peak_index + i
                    if min_ < self.cut_off_value:
                        break
                if v - x[min_index] > self.sensitivity_value:
                    break

            ub = min_index

        return lb, ub
