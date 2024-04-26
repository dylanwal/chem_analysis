from typing import Iterable

import numpy as np

from chem_analysis.processing.processing_method import ReSampling
from chem_analysis.processing.weigths.weights import Slices


class CutSlices(ReSampling):
    def __init__(self,
                 x_slices: slice | Iterable[slice] = None,  # TODO: generalize to n dimensions
                 y_slices: slice | Iterable[slice] = None,
                 invert: bool = False,
                 non_temporal_processing: bool = False
                 ):
        super().__init__(non_temporal_processing)
        if x_slices is None and y_slices is None:
            raise ValueError(f"Both '{type(self).__name__}.x_step' and '{type(self).__name__}.y_step' can't be None.")
        self.x_slices = x_slices
        self.y_slices = y_slices
        self.invert = invert

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.x_slices is None:
            raise ValueError(f"'{type(self).__name__}.x_slices' needs to be defined.")
        slice_ = Slices(self.x_slices, invert=self.invert)
        mask = slice_.get_mask(x, y)
        return x[mask], y[mask]

    def run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.x_slices is not None:
            slice_x = Slices(self.x_slices, invert=self.invert)
            mask_x = slice_x.get_mask(x, z)
        else:
            mask_x = np.ones_like(x, dtype=np.bool)

        if self.y_slices is not None:
            slice_y = Slices(self.y_slices, invert=self.invert)
            mask_y = slice_y.get_mask(y, z)
        else:
            mask_y = np.ones_like(y, dtype=np.bool)

        z = z[mask_y]
        return x[mask_x], y[mask_y], z[:, mask_x]

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError("this should never be called as 'run2D' is overloaded")

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
