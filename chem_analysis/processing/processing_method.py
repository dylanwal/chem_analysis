import abc
from typing import Any

import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.base_obj.signal_3d import Signal3D
from chem_analysis.utils.code_for_subclassing import MixinSubClassList


def get_class_instance_attributes(class_) -> list[tuple[str, Any]]:
    attrs = class_.__dict__.items()
    return [(k, v) for k, v in list(attrs) if v is not None]


def try_to_convert_to_number(text: str) -> int | float | str:
    try:
        num = float(text)
        if int(num) == float(num):
            return int(num)
        return num
    except ValueError:
        return text


class ProcessingMethod(MixinSubClassList, abc.ABC):

    def __init__(self, temporal_processing: int = 1):
        """

        Parameters
        ----------
        temporal_processing:
            0: data is not a time series
            1: data is a time series in 1 dimension
            2: data is a time series in 2 dimension
        """
        self.temporal_processing = temporal_processing

    def processing_str(self) -> str:
        """ generates string with method and args to record history of processing """
        text = type(self).__name__
        attrs = get_class_instance_attributes(self)
        text += "(" + ", ".join(f"{k}: {v}" for k, v in attrs) + ")"
        return text

    def parse_processing_str(self, text: str):  # -> ProcessingMethod
        """ experimental feature """
        method, args = text.split("(", maxsplit=1)
        sub_classes = self.all_sub_classes()
        for sub_class in sub_classes:
            if sub_class.__name__ == method:
                method = sub_class
                break
        else:
            raise ValueError(f"Method {method} not found.")

        args = args.replace(")", "").split(", ")
        kwargs = dict()
        for arg in args:
            try:
                k, v = arg.split(":", maxsplit=1)
                kwargs[k] = try_to_convert_to_number(v)
            except ValueError:
                raise ValueError(f"Argument {arg} is not a parsed successfully.")

        return method(**kwargs)

    def run(self, signal: Signal) -> Signal:
        x, y = self.run_xy(signal.x, signal.y)
        sig = signal.copy_with(x, y)
        sig.process_history.append(self.processing_str())
        return sig

    @abc.abstractmethod
    def run_xy(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ...

    def run2D(self, signal: Signal2D) -> Signal2D:
        x, y, z = self.run_xyz(signal.x, signal.y, signal.z)
        sig = signal.copy_with(x, y, z)
        sig.process_history.append(self.processing_str())
        return sig

    def run_xyz(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.temporal_processing >= 1:
            return self._run2D_temporal(x, y, z)
        return self._run2D(x, y, z)

    @abc.abstractmethod
    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ...

    def _run2D_temporal(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ y-axis is taken to be time series """
        for i in range(z.shape[0]):
            _, z[i, :] = self.run_xy(x, z[i, :])
        return x, y, z

    def run3D(self, signal: Signal3D) -> Signal3D:
        x, y, z, w = self.run_xyzw(signal.x, signal.y, signal.z, signal.w)
        sig = signal.copy_with(x, y, z, w)
        sig.process_history.append(self.processing_str())
        return sig

    def run_xyzw(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, w: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        if self.temporal_processing >= 1:
            self.temporal_processing -= 1  # reduce number of dimensions left as timeseries
            x, y, z, w = self._run3D_temporal(x, y, z, w)
            self.temporal_processing += 1  # reset original dimensions that are timeseries
            return x, y, z, w
        return self._run3D(x, y, z, w)

    @abc.abstractmethod
    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, w: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        ...
    
    def _run3D_temporal(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, w: np.ndarray) \
            -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """ z-axis is taken to be time series """
        for i in range(z.shape[0]):
            _, _, z[i, :] = self.run_xyz(x, y, z[i, :])
        return x, y, z, w


class Translation(ProcessingMethod, abc.ABC):
    ...


class Smoothing(ProcessingMethod, abc.ABC):
    ...

    def __call__(self, y: np.ndarray) -> np.ndarray:
        x, y = self.run_xy(np.arange(len(y)), y)
        return y


class Resampling(ProcessingMethod, abc.ABC):
    ...


class PhaseCorrection(ProcessingMethod, abc.ABC):
    ...


class FourierTransform(ProcessingMethod, abc.ABC):
    ...


class Edit(ProcessingMethod, abc.ABC):
    ...


class Baseline(ProcessingMethod, abc.ABC):
    def __init__(self,
                 temporal_processing: int = 1,
                 save_result: bool = False
                 ):
        super().__init__(temporal_processing)

        # for saving intermediate results
        self.save_result = save_result
        self.baseline = None
        self.x = None
        self.data = None

    def run_xy(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        baseline = self.get_baseline(x, y)
        data = y - baseline
        # data[y==0] = 0

        if self.save_result:
            self.baseline = baseline
            self.x = x
            self.data = y

        return x, data

    @abc.abstractmethod
    def get_baseline(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        ...

    def _run2D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()

    def _run3D(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError()
