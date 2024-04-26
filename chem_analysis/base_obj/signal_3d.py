from typing import Sequence, Iterable

import numpy as np

from chem_analysis.processing.base import Processor
from chem_analysis.analysis.peak import PeakBounded
from chem_analysis.base_obj.signal_2d import Signal2D


def validate_input(x_raw: np.ndarray, y_raw: np.ndarray, z_raw: np.ndarray, data_raw: np.ndarray):
    if len(x_raw.shape) != 1:
        raise ValueError(f"'x_raw' must shape 1. \n\treceived: {x_raw.shape}")
    if len(y_raw.shape) != 1:
        raise ValueError(f"'y_raw' must shape 1. \n\treceived: {y_raw.shape}")
    if len(z_raw.shape) != 1:
        raise ValueError(f"'z_raw' must shape 1. \n\treceived: {z_raw.shape}")
    if len(data_raw.shape) != 3:
        raise ValueError(f"'data_raw' must shape 2. \n\treceived: {data_raw.shape}")
    if x_raw.shape != data_raw.shape[2]:
        raise ValueError(f"'x_raw' and 'data_raw[2]' must have same shape. \n\treceived: x_raw:{x_raw.shape} "
                         f"|| data_raw: {data_raw.shape[2]}")
    if y_raw.shape != data_raw.shape[1]:
        raise ValueError(f"'y_raw' and 'data_raw[1]' must have same shape. \n\treceived: y_raw:{y_raw.shape} "
                         f"|| data_raw: {data_raw.shape[1]}")
    if z_raw.shape != data_raw.shape[1]:
        raise ValueError(f"'z_raw' and 'data_raw[1]' must have same shape. \n\treceived: z_raw:{z_raw.shape} "
                         f"|| data_raw: {data_raw.shape[1]}")


class Signal3D:
    """ signal 3D

    A signal is any x-y-z-zz data.

    """
    __count = 0
    _peak_type = PeakBounded

    def __init__(self,
                 x_raw: np.ndarray,
                 y_raw: np.ndarray,
                 z_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 zz_label: str = None,
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
        z_raw: np.ndarray[k]
            raw z data, length k
        data_raw: np.ndarray[k,j,i]
            raw z data, shape k,j,i
        x_label: str
            x-axis label
        y_label: str
            y-axis label
        z_label: str
            z-axis label
        zz_label: str
            4-axis label
        name: str
            user defined name
        """
        validate_input(x_raw, y_raw, z_raw, data_raw)

        self.x_raw = x_raw
        self.y_raw = y_raw
        self.z_raw = z_raw
        self.data_raw = data_raw
        self.id_ = id_ or Signal3D.__count
        Signal3D.__count += 1
        self.name = name or f"signal3D_{self.id_}"
        self.x_label = x_label or "x_axis"
        self.y_label = y_label or "y_axis"
        self.z_label = z_label or "z_axis"
        self.zz_label = zz_label or "zz_axis"

        self.processor = Processor()
        self._x = None
        self._y = None
        self._z = None
        self._data = None

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label}"
        text += f" (pts: {len(self)})"
        return text

    def __len__(self) -> int:
        return len(self.x)

    def _process(self):
        self._x, self._y, self._data = self.processor.run(self.x_raw, self.y_raw, self.data)

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
    def z(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._z

    @property
    def data(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()

        return self._data

    @property
    def number_of_signals(self):
        return len(self.y_raw)

    def pop(self, index: int) -> Signal2D:
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

    def get_signal(self, z_index: int, processed: bool = False) -> Signal2D:
        if processed:
            sig = Signal2D(x_raw=self.x, y_raw=self.y, data_raw=self.data[z_index, :, :], x_label=self.x_label,
                         y_label=self.y_label, name=f"slice_{self.z_label}: {self.z[z_index]}", id_=z_index)
        else:
            sig = Signal2D(x_raw=self.x_raw, y_raw=self.y_raw, data_raw=self.data_raw[z_index, :, :], x_label=self.x_label,
                         y_label=self.y_label, name=f"slice_{self.z_label}: {self.z[z_index]}", id_=z_index)
            sig.processor = self.processor.get_copy()
        sig.z_value = self.y[z_index]
        return sig

    @classmethod
    def from_signals(cls, signals: Sequence[Signal2D], z: np.ndarray = None, zz_label: str = None):  # -> Signal3D
        """ Turn Sequence of Signal2Ds into a Signal2D"""
        if z and len(z.shape) != 1 and z.shape[0] == len(signals):
            raise ValueError("The number of signals must be the same as the number of y points.\n"
                             f"\tnumber of signals: {len(signals)}\n\tnumber of z points:{z.shape[0]}")

        x = signals[0].x
        y = signals[0].y
        x_label = signals[0].x_label
        y_label = signals[0].y_label
        z_label = signals[0].z_label

        z = z or np.empty(len(signals))
        data = np.empty((len(signals), len(x)), dtype=signals[0].y.dtype)
        for i, sig in enumerate(signals):
            if np.all(sig.x != x):
                raise ValueError(f"Signal {i} has a different x-axis than first signal.")
            data[i, :] = sig.y
            if hasattr(sig, "time_"):
                z[i] = sig.time_
            else:
                z[i] = i

        return cls(x_raw=x, y_raw=y, z_raw=z, data_raw=data, x_label=x_label, y_label=y_label,
                   z_label=z_label, zz_label=z_label)

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################

    # @classmethod
    # def from_file(cls, path: str | pathlib.Path):
    #     from chem_analysis.utils.feather_format import feather_to_numpy
    #     from chem_analysis.utils.math import unpack_signal2D
    #
    #     if isinstance(path, str):
    #         path = pathlib.Path(path)
    #
    #     elif path.suffix == ".npy":
    #         data = np.load(str(path))
    #         x, y, z, data = unpack_signal2D(data)
    #         x_label = y_label = z_label = None
    #     else:
    #         raise NotImplemented("File type currently not supported.")
    #
    #     return cls(x_raw=x, y_raw=y, data_raw=data, x_label=x_label, y_label=y_label, z_label=z_label)
    #
    # def to_feather(self, path: str | pathlib.Path):
    #     from chem_analysis.utils.feather_format import numpy_to_feather
    #     from chem_analysis.utils.math import pack_time_series
    #
    #     headers = list(str(0) for i in range(len(self.y)+1))
    #     headers[0] = self.x_label
    #     headers[1] = self.y_label
    #     headers[2] = self.z_label
    #
    #     numpy_to_feather(pack_time_series(self.x, self.y, self.data), path, headers=headers)
    #
    # def to_csv(self, path: str | pathlib.Path, **kwargs):
    #     from chem_analysis.utils.math import pack_time_series
    #
    #     if "encodings" not in kwargs:
    #         kwargs["encoding"] = "utf-8"
    #     if "delimiter" not in kwargs:
    #         kwargs["delimiter"] = ","
    #
    #     np.savetxt(path, pack_time_series(self.x, self.time, self.data), **kwargs)  # noqa
    #
    # def to_npy(self, path: str | pathlib.Path, **kwargs):
    #     from chem_analysis.utils.math import pack_time_series
    #
    #     np.save(path, pack_time_series(self.x, self.y, self.data), **kwargs)
