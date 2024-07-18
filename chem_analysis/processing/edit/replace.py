from typing import Iterable, Sequence

import numpy as np

from chem_analysis.processing.processing_method import Edit
from chem_analysis.processing.weigths.weights import Slices, Spans


class ReplaceSlices(Edit):
    def __init__(self,
                 value: int | float = None,
                 x_slices: slice | Iterable[slice] = None,  # TODO: generalize to n dimensions
                 y_slices: slice | Iterable[slice] = None,
                 invert: bool = False,
                 temporal_processing: int = 1
                 ):
        super().__init__(temporal_processing)
        if x_slices is None and y_slices is None:
            raise ValueError(f"Both '{type(self).__name__}.x_step' and '{type(self).__name__}.y_step' can't be None.")
        self.x_slices = x_slices
        self.y_slices = y_slices
        self.invert = invert
        self.value = value

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.x_slices is None:
            raise ValueError(f"'{type(self).__name__}.x_slices' needs to be defined.")
        slice_ = Slices(self.x_slices, invert=self.invert)
        mask = slice_.get_mask(x, y)
        y[mask] = self.value
        return x, y

    def run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.x_slices is not None:
            slice_x = Slices(self.x_slices, invert=self.invert)
            mask_x = slice_x.get_mask(x, z)
        else:
            mask_x = np.ones_like(x, dtype=np.bool_)

        if self.y_slices is not None:
            slice_y = Slices(self.y_slices, invert=self.invert)
            mask_y = slice_y.get_mask(y, z)
        else:
            mask_y = np.ones_like(y, dtype=np.bool_)

        z[mask_y, mask_x] = self.value
        return x, y, z

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError("this should never be called as 'run2D' is overloaded")

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()


class ReplaceSpans(Edit):
    def __init__(self,
                 value: int | float,
                 x_spans: Sequence[float] | Iterable[Sequence[float]] = None,  # Sequence of length 2  # TODO: generalize to n dimensions
                 y_spans: Sequence[float] | Iterable[Sequence[float]] = None,  # Sequence of length 2
                 invert: bool = False,
                 temporal_processing: int = 1
                 ):
        super().__init__(temporal_processing)
        if x_spans is None and y_spans is None:
            raise ValueError("Both 'EveryN.x_step' and 'EveryN.y_step' can't be None.")
        self.x_spans = x_spans
        self.y_spans = y_spans
        self.invert = invert
        self.value = value

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.x_spans is None:
            raise ValueError(f"'{type(self).__name__}.x_spans' needs to be defined.")
        slice_ = Spans(self.x_spans, invert=self.invert)
        mask = slice_.get_mask(x, y)
        y[mask] = self.value
        return x, y

    def run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.x_spans is not None:
            slice_x = Spans(self.x_spans, invert=self.invert)
            mask_x = slice_x.get_mask(x, z)
            z[:, mask_x] = self.value

        if self.y_spans is not None:
            slice_y = Spans(self.y_spans, invert=self.invert)
            mask_y = slice_y.get_mask(y, z)
            z[mask_y, :] = self.value

        return x, y, z

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError("this should never be called as 'run2D' is overloaded")

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

