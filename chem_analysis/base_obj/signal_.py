from typing import Sequence

import numpy as np

import chem_analysis.algorithms.general_math as general_math
from chem_analysis.algorithms.processing.base import Processor


class Signal:
    """ signal

    A signal is any x-y data.

    Attributes
    ----------
    name: str
        Any name the user wants to add.
    x_label: str
        x-axis label
    y_label: str
        y-axis label
    """
    __count = 0

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 name: str = None,
                 id_: int = None
                 ):
        """

        Parameters
        ----------
        x: np.ndarray
            x data
        y: np.ndarray
            y data
        x_label: str
            x-axis label
        y_label: str
            y-axis label
        name: str
            user defined name

        Notes
        -----
        * Either 'ser' or 'x' and 'y' are required but not both.

        """
        self.x_raw = x
        self.y_raw = y
        self.id_ = id_ if id_ is not None else Signal.__count
        Signal.__count += 1
        self.name = name if name is not None else f"signal_{self.id_}"
        self.x_label = x_label if x_label is not None else "x_axis"
        self.y_label = y_label if y_label is not None else "y_axis"

        self.processor = Processor()
        self._x = None
        self._y = None

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self)})"
        return text

    def __len__(self) -> int:
        return len(self.x)

    @property
    def x(self) -> np.ndarray:
        if not self.processor.processed:
            self._x, self._y = self.processor.run(self.x_raw, self.y_raw)

        return self._x

    @property
    def y(self) -> np.ndarray:
        if not self.processor.processed:
            self._x, self._y = self.processor.run(self.x_raw, self.y_raw)

        return self._y

    @property
    def y_normalized_by_max(self, x_range: Sequence[int | float] = None) -> np.ndarray:
        if x_range is None:
            return self.y/np.max(self.y)

        slice_ = general_math.get_slice(self.x, *x_range)
        return general_math.normalize_by_max(self.y[slice_])

    @property
    def y_normalized_by_area(self, x_range: Sequence[int | float] = None) -> np.ndarray:
        if x_range is None:
            return general_math.normalize_by_area(self.x, self.y)

        slice_ = general_math.get_slice(self.x, *x_range)
        return general_math.normalize_by_area(self.x[slice_], self.y[slice_])
