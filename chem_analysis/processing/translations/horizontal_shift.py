from typing import Sequence

import numpy as np

from chem_analysis.processing.processing_method import Translation


class HorizontalShift(Translation):
    def __init__(self,
                 shift_x: float | Sequence[float] = None,
                 shift_index: int | Sequence[int] = None,
                 wrap: bool = True,
                 temporal_processing: int = 1
                 ):
        super().__init__(temporal_processing)
        if (shift_index is None) == (shift_x is None):
            raise ValueError("Provide only one: shift_x or shift_index")
        self.shift_x = shift_x
        self.shift_index = shift_index
        self.wrap = wrap

    def _get_shift_index(self, x: np.ndarray):
        if self.shift_index is not None:
            return

        if isinstance(self.shift_x, Sequence):
            self.shift_index = np.zeros_like(self.shift_x, dtype=np.int32)
            for i in range(len(self.shift_x)):
                self.shift_index[i] = self._get_shift_index_single(x, self.shift_x[i])
        else:
            self.shift_index = self._get_shift_index_single(x, self.shift_x)

    @staticmethod
    def _get_shift_index_single(x: np.ndarray, shift_x: float) -> int:
        if abs(shift_x) > (np.max(x) - np.min(x)):
            raise ValueError("shift is larger than x range")
        if shift_x > 0:
            new_min = np.min(x) + shift_x
            index = np.argmin(np.abs(x - new_min))
        else:
            new_max = np.max(x) + shift_x
            index = np.argmin(np.abs(x - new_max)) - len(x)
        return index

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if isinstance(self.shift_index, Sequence) or isinstance(self.shift_x, Sequence):
            raise ValueError("For x-y signals, provide single values for 'shift_index' and 'shift_x'.")
        self._get_shift_index(x)

        if self.shift_index == 0:
            return x, y

        if self.wrap:
            y = np.roll(y, self.shift_index)
            return x, y

        if self.shift_index > 0:
            y = y[:-self.shift_index]
            x = x[self.shift_index:]
        else:
            y = y[self.shift_index:]
            x = x[:-self.shift_index]
        return x, y

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self._get_shift_index(x)

        if isinstance(self.shift_index, int):
            if self.shift_index == 0:
                return x, y, z
            if self.wrap:
                z = np.roll(z, self.shift_index, axis=1)
                return x, y, z
            if self.shift_index > 0:
                z = z[:, -self.shift_index]
                x = x[self.shift_index:]
            else:
                z = z[:, self.shift_index:]
                x = x[:-self.shift_index]
            return x, y, z

        if not self.wrap:
            raise ValueError("'wrap must be true otherwise different x-axis are needed.")
        for i in range(z.shape[0]):
            z[i, :] = np.roll(z[i, :], self.shift_index[i])

        return x, y, z

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
