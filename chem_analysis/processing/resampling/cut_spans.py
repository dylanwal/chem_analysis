from typing import Sequence, Iterable

import numpy as np

from chem_analysis.processing.processing_method import ReSampling
from chem_analysis.processing.weigths.weights import Spans


class CutSpans(ReSampling):
    def __init__(self,
                 x_spans: Sequence[float] | Iterable[Sequence[float]] = None,  # Sequence of length 2  # TODO: generalize to n dimensions
                 y_spans: Sequence[float] | Iterable[Sequence[float]] = None,  # Sequence of length 2
                 invert: bool = False,
                 non_temporal_processing: bool = False
                 ):
        super().__init__(non_temporal_processing)
        if x_spans is None and y_spans is None:
            raise ValueError("Both 'EveryN.x_step' and 'EveryN.y_step' can't be None.")
        self.x_spans = x_spans
        self.y_spans = y_spans
        self.invert = invert

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.x_spans is None:
            raise ValueError(f"'{type(self).__name__}.x_spans' needs to be defined.")
        slice_ = Spans(self.x_spans, invert=self.invert)
        mask = slice_.get_mask(x, y)
        return x[mask], y[mask]

    def run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.x_spans is not None:
            slice_x = Spans(self.x_spans, invert=self.invert)
            mask_x = slice_x.get_mask(x, z)
        else:
            mask_x = np.ones_like(x, dtype=np.bool)

        if self.y_spans is not None:
            slice_y = Spans(self.y_spans, invert=self.invert)
            mask_y = slice_y.get_mask(y, z)
        else:
            mask_y = np.ones_like(y, dtype=np.bool)

        return x[mask_x], y[mask_y], z[mask_y, mask_x]

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError("this should never be called as 'run2D' is overloaded")

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()