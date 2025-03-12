from typing import Sequence

import numpy as np

from chem_analysis.processing.processing_method import ProcessingMethod, Baseline


class CompoundProcessingBaseline(Baseline):
    def __init__(self,
                 methods: Sequence[ProcessingMethod],
                 save_result: bool = False,
                 temporal_processing: int = 1,
                 ):
        """
        Compound Processing Baseline Class
        allows for smoothing, masking, subsampling, etc. prior to computing baseline
        the methods will not affect x,y signal; just the data used for computing baseline

        Parameters
        ----------
        methods:
            list of processing methods applied to x,y before computing baseline
            last method in list must be a 'Baseline' ProcessingMethod
        temporal_processing:

        """
        if not isinstance(methods[-1], Baseline):
            raise ValueError("The last method must be a 'Baseline' processing")

        super().__init__(temporal_processing, save_result)
        self.methods = methods

    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        x_, y_ = np.copy(x), np.copy(y)
        for method in self.methods[:-1]:
           x_, y_ = method.run(x, y)
        baseline = self.methods[-1].get_baseline(x_,y_)
        return np.interp(x, x_, baseline)

    # def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    #
    #
    #     if self.save_result:
    #         self.baseline = baseline
    #         self.x = x
    #         self.data = y
    #
    #
    #     return x, y - baseline
    #
    # def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    #     pass
    #
    # def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, w: np.ndarray) -> tuple[
    #     np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    #     pass