import pathlib

import numpy as np

from chem_analysis.processing.base import Processor
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.peak import PeakBoundedStats


class SignalArray:
    """
    A grouping of Signals where each Signal occurred at a different time interval.
    """
    _peak_type = PeakBoundedStats

    def __init__(self,
                 x_raw: np.ndarray,
                 time_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None
                 ):
        self.name = name
        self.x_raw = x_raw
        self.time_raw = time_raw
        self.raw_data = data_raw
        self.x_label = x_label if x_label is not None else "x_axis"
        self.y_label = y_label if y_label is not None else "time"
        self.z_label = z_label if z_label is not None else "z_axis"

        self.processor = Processor()
        self._x = None
        self._time = None
        self._data = None

    def _process(self):
        self._x, self._time, self._data = self.processor.run(self.x_raw, self.time_raw, self.raw_data)

    @property
    def x(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._x

    @property
    def time(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._time

    @property
    def time_zeroed(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._time - self.time_raw[0]

    @property
    def data(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._data

    @property
    def number_of_signals(self):
        return len(self.time_raw)

    def get_signal(self, index: int) -> Signal:
        sig = Signal(x_raw=self.x, y_raw=self.data[index, :], x_label=self.x_label,
                     y_label=self.y_label, name=f"time: {self.time[index]}", id_=index)
        sig.processor = self.processor.get_copy()
        return sig

    @classmethod
    def from_file(cls, path: str | pathlib.Path):
        from chem_analysis.utils.feather_format import feather_to_numpy, unpack_time_series
        if isinstance(path, str):
            path = pathlib.Path(path)

        if path.suffix == ".csv":
            data = np.loadtxt(path, delimiter=",")
            x, time_, data = unpack_time_series(data)
        elif path.suffix == ".feather":
            data = feather_to_numpy(path)
            x, time_, data = unpack_time_series(data)
        else:
            raise NotImplemented("File type currently not supported.")

        return cls(x, time_, data)
