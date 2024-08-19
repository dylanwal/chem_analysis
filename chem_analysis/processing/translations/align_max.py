
import numpy as np

from chem_analysis.processing.processing_method import Translation
from chem_analysis.utils.math import get_slice
from chem_analysis.processing.translations.horizontal_shift import HorizontalShift


class AlignMax(Translation):
    def __init__(self,
                 range_: tuple[float, float] = None,
                 range_index: slice = None,
                 x_value: int | float = None,
                 wrap: bool = True,
                 temporal_processing: int = 1
                 ):
        super().__init__(temporal_processing)
        self.range_ = range_
        self.range_index = range_index
        self.x_value = x_value
        self.wrap = wrap
        self.x_value_index = None
        self.shift_index = None

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.x_value is None:
            raise ValueError("Set the 'AlignMax.x_value' you want to align the max to.")
        if self.range_ is not None:
            self.range_index = get_slice(x, self.range_[0], self.range_[1])
        if self.range_index is None:
            self.range_index = slice(0, -1)

        # get shift index
        max_indices = np.argmax(y[self.range_index]) + self.range_index.start
        self.x_value_index = np.argmin(np.abs(x - self.x_value))
        self.shift_index = self.x_value_index - max_indices

        translation = HorizontalShift(shift_index=self.shift_index, wrap=self.wrap)
        return translation.run(x, y)

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.range_ is not None:
            self.range_index = get_slice(x, self.range_[0], self.range_[1])
        if self.range_index is None:
            self.range_index = slice(0, -1)
        if self.x_value is None:
            # use first spectra max as reference
            self.x_value_index = np.argmax(z[0, self.range_index]) + self.range_index.start
            self.x_value = x[self.x_value_index]
        else:
            self.x_value_index = np.argmin(np.abs(x - self.x_value))

        # get shift index
        max_indices = np.argmax(z[:, self.range_index], axis=1) + self.range_index.start
        self.shift_index = self.x_value_index - max_indices

        translation = HorizontalShift(shift_index=self.shift_index, wrap=self.wrap)
        return translation._run2D(x, y, z)

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
