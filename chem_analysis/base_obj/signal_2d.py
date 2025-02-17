import pathlib
from typing import Sequence, Iterable, Iterator

import numpy as np

from chem_analysis.base_obj.unify_methods import UnifyMethod, UnifyMethodStrict
from chem_analysis.processing.processor import Processor
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.utils.math import unpack_signal2D


def validate_input(x_raw: np.ndarray, y_raw: np.ndarray, z_raw: np.ndarray):
    if len(x_raw.shape) != 1:
        raise ValueError(f"'x_raw' must shape 1. \n\treceived: {x_raw.shape}")
    if len(y_raw.shape) != 1:
        raise ValueError(f"'y_raw' must shape 1. \n\treceived: {y_raw.shape}")
    if len(z_raw.shape) != 2:
        raise ValueError(f"'z_raw' must shape 2. \n\treceived: {z_raw.shape}")
    if x_raw.shape[0] != z_raw.shape[1]:
        raise ValueError(f"'x_raw' and 'z_raw[1]' must have same shape. \n\treceived: x_raw:{x_raw.shape} "
                         f"|| z_raw.shape[1]: {z_raw.shape[1]}")
    if y_raw.shape[0] != z_raw.shape[0]:
        raise ValueError(f"'y_raw' and 'z_raw[0]' must have same shape. \n\treceived: y_raw:{y_raw.shape} "
                         f"|| z_raw.shape[0]: {z_raw.shape[0]}")


class Signal2D:
    """ signal 2D

    A signal is any x-y-z data.

    """
    __count = 0
    _signal = Signal

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None,
                 id_: int = None
                 ):
        """

        Parameters
        ----------
        x: np.ndarray[i]
            raw x data, length i
        y: np.ndarray[j]
            raw y data, length j
        z: np.ndarray[j,i]
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
        validate_input(x, y, z)

        self.x_raw = x
        self.y_raw = y
        self.z_raw = z
        self.id_ = id_ or Signal2D.__count
        Signal2D.__count += 1
        self.name = name or f"signal_{self.id_}"
        self.x_label = x_label or "x_axis"
        self.y_label = y_label or "y_axis"
        self.z_label = z_label or "z_axis"

        self.processor = Processor()
        self._x = None
        self._y = None
        self._z = None

        self.extract_value = None  # value

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs. {self.y_label} vs. {self.z_label}"
        text += f" (shape: {self.z_raw.shape})"
        return text

    def _process(self):
        self._x, self._y, self._z = self.processor.run(self.x_raw, self.y_raw, self.z_raw)

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
    def number_of_signals(self):
        return len(self.y_raw)

    def pop(self, index: int) -> Signal:
        sig = self.get_signal(index)
        self.delete(index)
        return sig

    def delete(self, index: int | Iterable[int] | slice):
        if isinstance(index, int):
            index = [index]
        index.sort(reverse=True)  # delete largest to smallest to avoid issue of changing index
        for i in index:
            self.z_raw = np.delete(self.z_raw, i, axis=0)
            self.y_raw = np.delete(self.y_raw, i)

    def get_signal(self, y_index: int, processed: bool = True, copy_: bool = False) -> Signal:
        """

        Parameters
        ----------
        y_index
        processed:
            True: get x, z
            False: get x_raw, z_raw
        copy_:
            True: data will be a copy.
            False: data will be a view (until edited)

        Returns
        -------

        Should return a 'view' and not 'copy'. But will become a copy if edited.
        https://numpy.org/doc/stable/user/basics.copies.html

        """
        if processed:
            x, y = self.x, self.z[y_index, :]
        else:
            x, y = self.x_raw, self.z_raw[y_index, :]

        if copy_:
            x, y = x.copy(), y.copy()

        sig = self._signal(x=x, y=y, x_label=self.x_label, y_label=self.y_label,
                           name=f"slice_{self.y_label}: {self.y[y_index]}", id_=y_index)
        sig.extract_value = self.y[y_index]
        if not processed:
            sig.processor = self.processor.get_copy()
            sig.processor = False
        return sig

    def signal_iter(self) -> Iterator[Signal]:
        for i in range(self.number_of_signals):
            yield self.get_signal(i, processed=True)

    @classmethod
    def from_signals(cls,
                     signals: Sequence[Signal],
                     y: np.ndarray = None,
                     y_label: str = None,
                     unify_method: UnifyMethod = UnifyMethodStrict()
                     ):  # -> Signal2D
        """ Turn Sequence of Signals into a Signal2D"""
        if y is None:
            y = np.arange(len(signals))
        else:
            if len(y.shape) != 1 and y.shape[0] == len(signals):
                raise ValueError("The number of signals must be the same as the number of y points.\n"
                                 f"\tnumber of signals: {len(signals)}\n\tnumber of y points:{y.shape[0]}")

        x_label = signals[0].x_label
        z_label = signals[0].y_label
        x, z = unify_method.run(signals)
        return cls(x, y, z, x_label=x_label, y_label=y_label, z_label=z_label)

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################
    def to_feather(self, path: str | pathlib.Path):
        from chem_analysis.utils.feather_format import numpy_to_feather
        from chem_analysis.utils.math import pack_time_series

        headers = list(str(0) for i in range(len(self.y) + 1))
        headers[0] = self.x_label
        headers[1] = self.y_label
        headers[2] = self.z_label

        numpy_to_feather(pack_time_series(self.x, self.y, self.z), path, headers=headers)

    def to_csv(self, path: str | pathlib.Path, **kwargs):
        from chem_analysis.utils.math import pack_time_series

        if "encodings" not in kwargs:
            kwargs["encoding"] = "utf-8"
        if "delimiter" not in kwargs:
            kwargs["delimiter"] = ","

        np.savetxt(path, pack_time_series(self.x, self.time, self.z), **kwargs)  # noqa

    def to_npy(self, path: str | pathlib.Path, **kwargs):
        """Save an array to a binary file in NumPy ``.npy`` format."""
        from chem_analysis.utils.math import pack_time_series

        np.save(path, pack_time_series(self.x, self.y, self.z), **kwargs)

    def to_npz(self, path: str | pathlib.Path, **kwargs):
        """Save an array to a binary file in NumPy ``.npz`` format."""
        np.savez(path, x=self.x, y=self.y, z=self.z, **kwargs)

    @classmethod
    def from_file(cls, path: str | pathlib.Path):
        if isinstance(path, str):
            path = pathlib.Path(path)
        if path.suffix == ".npz":
            return cls.from_npz(path)
        elif path.suffix == ".feather":
            return cls.from_feather(path)
        elif path.suffix == ".csv":
            return cls.from_csv(path)
        elif path.suffix == ".npy":
            return cls.from_npy(path)
        else:
            raise ValueError(f"Unsupported file type: {path.suffix}")

    @classmethod
    def from_csv(cls, path: str | pathlib.Path):
        z = np.loadtxt(path, delimiter=",")
        x, y, z = unpack_signal2D(z)
        x_label = y_label = z_label = None
        return cls(x, y, z, x_label=x_label, y_label=y_label, z_label=z_label)

    @classmethod
    def from_feather(cls, path: str | pathlib.Path):
        from chem_analysis.utils.feather_format import feather_to_numpy
        z, names = feather_to_numpy(path)
        x, y, z = unpack_signal2D(z)
        if names[0] != "0":
            x_label = names[0]
            y_label = names[1]
            z_label = names[2]
        else:
            x_label = y_label = z_label = None
        return cls(x, y, z, x_label=x_label, y_label=y_label, z_label=z_label)

    @classmethod
    def from_npy(cls, path: str | pathlib.Path):
        z = np.load(str(path))
        x, y, z = unpack_signal2D(z)
        x_label = y_label = z_label = None
        return cls(x, y, z, x_label=x_label, y_label=y_label, z_label=z_label)

    @classmethod
    def from_npz(cls, path: str | pathlib.Path):
        npzfile = np.load(str(path))
        x, y, z = npzfile['x'], npzfile['y'], npzfile['z']
        x_label = y_label = z_label = None
        return cls(x, y, z, x_label=x_label, y_label=y_label, z_label=z_label)
