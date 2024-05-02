
import numpy as np

from chem_analysis.processing.processing_method import Translation
from chem_analysis.utils.math import get_slice


class ScaleMax(Translation):
    def __init__(self,
                 range_: tuple[float, float] = None,
                 range_index: slice = None,
                 new_max_value: int | float = None,
                 wrap: bool = True,
                 non_temporal_processing: bool = False
                 ):
        super().__init__(non_temporal_processing)
        self.range_ = range_
        self.range_index = range_index
        self.new_max_value = new_max_value
        self.wrap = wrap
        self.scale = None

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.new_max_value is None:
            raise ValueError("Set the 'AlignMax.x_value' you want to align the max to.")
        if self.range_ is not None:
            self.range_index = get_slice(x, self.range_[0], self.range_[1])
        if self.range_index is None:
            self.range_index = slice(0, -1)

        self.scale = self.new_max_value / np.max(y[self.range_index])
        return x, y * self.scale

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.range_ is not None:
            self.range_index = get_slice(x, self.range_[0], self.range_[1])
        if self.range_index is None:
            self.range_index = slice(0, -1)
        if self.new_max_value is None:
            # use first spectra max as reference
            self.new_max_value = np.max(z[0, self.range_index])

        # get shift index
        self.scale = self.new_max_value / np.max(z[:, self.range_index], axis=1)
        return x, y, z * self.scale.reshape(-1, 1)

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
