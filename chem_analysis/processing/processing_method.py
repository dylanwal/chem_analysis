import abc

import numpy as np

from chem_analysis.utils.code_for_subclassing import MixinSubClassList


class ProcessingMethod(MixinSubClassList, abc.ABC):

    def __init__(self, non_temporal_processing: bool = False):
        self.non_temporal_processing = non_temporal_processing

    @abc.abstractmethod
    def run(self, x: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...

    def run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.non_temporal_processing:
            return self._run2D(x, y, data)
        return self._run2D_temporal(x, y, data)

    def _run2D(self, x: np.ndarray, y: np.ndarray, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def _run2D_temporal(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        for i in range(z.shape[0]):
            _, z[i, :] = self.run(x, z[i, :])
        return x, y, z

    def run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()


class Translation(ProcessingMethod, abc.ABC):
    ...


class Smoothing(ProcessingMethod, abc.ABC):
    ...


class ReSampling(ProcessingMethod, abc.ABC):
    ...


class PhaseCorrection(ProcessingMethod, abc.ABC):
    ...


class FourierTransform(ProcessingMethod, abc.ABC):
    ...


class Edit(ProcessingMethod, abc.ABC):
    ...


class Baseline(ProcessingMethod, abc.ABC):
    def __init__(self,
                 non_temporal_processing: bool = False,
                 save_result: bool = False
                 ):
        super().__init__(non_temporal_processing)

        # for saving intermediate results
        self.save_result = save_result
        self.baseline = None
        self.x = None
        self.data = None

    def run(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        baseline = self.get_baseline(x, y)
        data = y - baseline

        if self.save_result:
            self.baseline = baseline
            self.x = x
            self.data = y

        return x, data

    @abc.abstractmethod
    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        ...
