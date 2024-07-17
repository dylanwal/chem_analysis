from typing import Sequence

import numpy as np

from chem_analysis.processing.processing_method import Resampling


class EveryN(Resampling):
    def __init__(self,
                 step: Sequence[int],
                 start_index: Sequence[int] = None,
                 temporal_processing: int = 1
                 ):
        """

        Parameters
        ----------
        step:
            the steps to take in resampling
            each position in the sequence corresponds to a dimension
        start_index:
            starting index
            each position in the sequence corresponds to a dimension
        """
        super().__init__(temporal_processing)
        if not any(isinstance(i, int) and i >= 0 for i in step):
            raise ValueError(f"'{type(self).__name__}.step' must be be positive(>=0) integers.")
        if start_index is not None:
            if not any(isinstance(i, int) and i >= 0 for i in start_index):
                raise ValueError(f"'{type(self).__name__}.start_index' must be be positive(>=0) integers.")
            if len(step) != len(start_index):
                raise ValueError(f"'{type(self).__name__}.step' and '{type(self).__name__}.start_index' must "
                                 f"have the same shape.")
        self.step = step
        self.start_index = start_index

    def run(self, x: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if len(self.step) != 1:
            raise ValueError(f"'{type(self).__name__}.step' needs to length 1.")
        if self.start_index is None:
            self.start_index = [0]
        else:
            if len(self.start_index) != 1:
                raise ValueError(f"'{type(self).__name__}.start_index' needs to length 1.")
            if self.start_index[0] > len(x):
                raise ValueError(f"'{type(self).__name__}.start_index' is larger then the length of the array")

        return x[self.start_index[0]::self.step[0]], data[self.start_index[0]::self.step[0]]

    def run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if len(self.step) != 2:
            raise ValueError(f"'{type(self).__name__}.step' needs to length 2.")
        if self.start_index is None:
            self.start_index = [0, 0]
        else:
            if len(self.start_index) != 2:
                raise ValueError(f"'{type(self).__name__}.start_index' needs to length 2.")
            if any((
                    self.start_index[0] > len(x),
                    self.start_index[1] > len(y)
            )):
                raise ValueError(f"'{type(self).__name__}.start_index' is larger then the length of the array")

        return x[self.start_index[0]::self.step[0]], y[self.start_index[1]::self.step[1]], \
               data[self.start_index[1]::self.step[1], self.start_index[0]::self.step[0]]

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError("this should never be called as 'run2D' is overloaded")

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        if len(self.step) != 3:
            raise ValueError(f"'{type(self).__name__}.step' needs to length 3.")
        if self.start_index is None:
            self.start_index = [0, 0, 0]
        else:
            if len(self.start_index) != 3:
                raise ValueError(f"'{type(self).__name__}.start_index' needs to length 3.")
            if any((
                    self.start_index[0] > len(x),
                    self.start_index[1] > len(y),
                    self.start_index[2] > len(z)
            )):
                raise ValueError(f"'{type(self).__name__}.start_index' is larger then the length of the array")

        return x[self.start_index[0]::self.step[0]], y[self.start_index[1]::self.step[1]], \
               y[self.start_index[2]::self.step[2]], \
               data[
               self.start_index[2]::self.step[2],
               self.start_index[1]::self.step[1],
               self.start_index[0]::self.step[0]
               ]
