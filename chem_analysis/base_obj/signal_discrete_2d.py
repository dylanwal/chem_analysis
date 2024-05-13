import pathlib
from typing import Sequence, Iterable

import numpy as np

from chem_analysis.processing.processor import Processor
from chem_analysis.analysis.peak import PeakDiscrete
from chem_analysis.base_obj.signal_discrete import SignalDiscrete


def validate_input(x_raw: np.ndarray, y_raw: np.ndarray, data_raw: np.ndarray):
    if len(x_raw.shape) != 1:
        raise ValueError(f"'x_raw' must shape 1. \n\treceived: {x_raw.shape}")
    if len(y_raw.shape) != 1:
        raise ValueError(f"'y_raw' must shape 1. \n\treceived: {y_raw.shape}")
    if len(data_raw.shape) != 2:
        raise ValueError(f"'data_raw' must shape 2. \n\treceived: {data_raw.shape}")
    if x_raw.shape[0] != data_raw.shape[1]:
        raise ValueError(f"'x_raw' and 'data_raw[1]' must have same shape. \n\treceived: x_raw:{x_raw.shape} "
                         f"|| data_raw: {data_raw.shape[1]}")
    if y_raw.shape[0] != data_raw.shape[0]:
        raise ValueError(f"'y_raw' and 'data_raw[0]' must have same shape. \n\treceived: y_raw:{y_raw.shape} "
                         f"|| data_raw: {data_raw.shape[0]}")


class SignalDiscrete2D:
    """ signal 2D

    A signal is any x-y-z data.

    """
    __count = 0
    _peak_type = PeakDiscrete
    _signal = SignalDiscrete

    def __init__(self,
                 x_raw: np.ndarray,
                 y_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None,
                 id_: int = None
                 ):
        """

        Parameters
        ----------
        x_raw: np.ndarray[i]
            raw x data, length i
        y_raw: np.ndarray[j]
            raw y data, length j
        data_raw: np.ndarray[j,i]
            raw z data, shape j,i
        x_label: str
            x-axis label
        y_label: str
            y-axis label
        z_label: str
            z-axis label
        name: str
            user defined name
        """
        validate_input(x_raw, y_raw, data_raw)

        self.x_raw = x_raw
        self.y_raw = y_raw
        self.data_raw = data_raw
        self.id_ = id_ or SignalDiscrete2D.__count
        SignalDiscrete2D.__count += 1
        self.name = name or f"signal_{self.id_}"
        self.x_label = x_label or "x_axis"
        self.y_label = y_label or "y_axis"
        self.z_label = z_label or "z_axis"

        self.processor = Processor()
        self._x = None
        self._y = None
        self._data = None

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs. {self.y_label} vs. {self.z_label}"
        text += f" (pts: {len(self)})"
        return text

    def __len__(self) -> int:
        return len(self.y)

    def _process(self):
        self._x, self._y, self._data = self.processor.run(self.x_raw, self.y_raw, self.data_raw)

    @property
    def x(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._x

    @property
    def y(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._y

    @property
    def data(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._data

    @property
    def number_of_signals(self):
        return len(self.y_raw)

    def pop(self, index: int) -> SignalDiscrete:
        sig = self.get_signal(index)
        self.delete(index)
        return sig

    def delete(self, index: int | Iterable):
        if isinstance(index, int):
            index = [index]
        index.sort(reverse=True)  # delete largest to smallest to avoid issue of changing index
        for i in index:
            self.data_raw = np.delete(self.data_raw, i, axis=0)
            self.y_raw = np.delete(self.y_raw, i)

    def get_signal(self, y_index: int, processed: bool = False) -> SignalDiscrete:
        if processed:
            sig = self._signal(x_raw=self.x, data_raw=self.data[y_index, :], x_label=self.x_label,
                         y_label=self.y_label, name=f"slice_{self.y_label}: {self.y[y_index]}", id_=y_index)
        else:
            sig = self._signal(x_raw=self.x_raw, data_raw=self.data_raw[y_index, :], x_label=self.x_label,
                         y_label=self.y_label, name=f"slice_{self.y_label}: {self.y[y_index]}", id_=y_index)
            sig.processor = self.processor.get_copy()
        sig.y_value = self.y[y_index]
        return sig

    @classmethod
    def from_signals(cls, signals: Sequence[SignalDiscrete], y: np.ndarray = None):  # -> Signal2D
        """ Turn Sequence of Signals into a Signal2D"""
        # TODO: add interpolation option if x-axis not same
        if y is not None and len(y.shape) != 1 and y.shape[0] == len(signals):
            raise ValueError("The number of signals must be the same as the number of y points.\n"
                             f"\tnumber of signals: {len(signals)}\n\tnumber of y points:{y.shape[0]}")

        x = signals[0].x
        x_label = signals[0].x_label
        z_label = signals[0].y_label

        if y is None:
            y = np.empty(len(signals))
        data = np.empty((len(signals), len(x)), dtype=signals[0].y.dtype)
        for i, sig in enumerate(signals):
            if not np.all(np.isclose(sig.x, x, rtol=0.01)):
                raise ValueError(f"Signal {i} has a different x-axis than first signal.")
            data[i, :] = sig.y
            if hasattr(sig, "time_"):
                y[i] = sig.time_
            else:
                y[i] = i

        return cls(x_raw=x, y_raw=y, data_raw=data, x_label=x_label, z_label=z_label)

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################

    @classmethod
    def from_file(cls, path: str | pathlib.Path):
        from chem_analysis.utils.feather_format import feather_to_numpy
        from chem_analysis.utils.math import unpack_signal2D

        if isinstance(path, str):
            path = pathlib.Path(path)

        if path.suffix == ".csv":
            data = np.loadtxt(path, delimiter=",")
            x, y, data = unpack_signal2D(data)
            x_label = y_label = z_label = None

        elif path.suffix == ".feather":
            data, names = feather_to_numpy(path)
            x, y, data = unpack_signal2D(data)
            if names[0] != "0":
                x_label = names[0]
                y_label = names[1]
                z_label = names[2]
            else:
                x_label = y_label = z_label = None

        elif path.suffix == ".npy":
            data = np.load(str(path))
            x, y, data = unpack_signal2D(data)
            x_label = y_label = z_label = None
        else:
            raise NotImplemented("File type currently not supported.")

        return cls(x_raw=x, y_raw=y, data_raw=data, x_label=x_label, y_label=y_label, z_label=z_label)

    def to_feather(self, path: str | pathlib.Path):
        from chem_analysis.utils.feather_format import numpy_to_feather
        from chem_analysis.utils.math import pack_time_series

        headers = list(str(0) for i in range(len(self.y)+1))
        headers[0] = self.x_label
        headers[1] = self.y_label
        headers[2] = self.z_label

        numpy_to_feather(pack_time_series(self.x, self.y, self.data), path, headers=headers)

    def to_csv(self, path: str | pathlib.Path, **kwargs):
        from chem_analysis.utils.math import pack_time_series

        if "encodings" not in kwargs:
            kwargs["encoding"] = "utf-8"
        if "delimiter" not in kwargs:
            kwargs["delimiter"] = ","

        np.savetxt(path, pack_time_series(self.x, self.time, self.data), **kwargs)  # noqa

    def to_npy(self, path: str | pathlib.Path, **kwargs):
        """Save an array to a binary file in NumPy ``.npy`` format."""
        from chem_analysis.utils.math import pack_time_series

        np.save(path, pack_time_series(self.x, self.y, self.data), **kwargs)
