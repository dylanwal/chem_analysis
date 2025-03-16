import pathlib
from typing import Sequence, Iterable, Iterator
import copy

import numpy as np

from chem_analysis.base_obj.parameters import Parameters
from chem_analysis.base_obj.unify_methods import UnifyMethod, UnifyMethodStrict
from chem_analysis.base_obj.signal_ import Signal, check_array_inf_nan
from chem_analysis.utils.math import unpack_signal2D


def validate_input(x: np.ndarray, y: np.ndarray, z: np.ndarray):
    if len(x.shape) != 1:
        raise ValueError(f"'x' must shape 1. \n\treceived: {x.shape}")
    if len(y.shape) != 1:
        raise ValueError(f"'y' must shape 1. \n\treceived: {y.shape}")
    if len(z.shape) != 2:
        raise ValueError(f"'z' must shape 2. \n\treceived: {z.shape}")
    if x.shape[0] != z.shape[1]:
        raise ValueError(f"'x' and 'z[1]' must have same shape. \n\treceived: x:{x.shape} "
                         f"|| z.shape[1]: {z.shape[1]}")
    if y.shape[0] != z.shape[0]:
        raise ValueError(f"'y' and 'z[0]' must have same shape. \n\treceived: y:{y.shape} "
                         f"|| z.shape[0]: {z.shape[0]}")
    check_array_inf_nan(x, name="x")
    check_array_inf_nan(y, name="y")
    check_array_inf_nan(y, name="z")


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
                 id_: int = None,
                 parameters: Parameters = None,
                 process_history: bool = False,
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
        parameters: Parameters
            various meta-data
        process_history: list[str]
            List of processing methods the signal has been through
        """
        validate_input(x, y, z)

        self.x = x
        self.y = y
        self.z = z
        self.id_ = id_ or Signal2D.__count
        Signal2D.__count += 1
        self.name = name or f"signal_{self.id_}"
        self.x_label = x_label or "x_axis"
        self.y_label = y_label or "y_axis"
        self.z_label = z_label or "z_axis"

        self.parameters = parameters
        self.extract_value = None

        if isinstance(process_history, str):
            process_history = [process_history]
        self.process_history = process_history or []

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs. {self.y_label} vs. {self.z_label}"
        text += f" (shape: {self.z.shape})"
        return text

    def copy_with(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, deep: bool = True):
        if deep:
            copy_method = copy.deepcopy
        else:
            copy_method = copy.copy

        return Signal2D(
            x,
            y,
            z,
            name=copy_method(self.name),
            x_label=copy_method(self.x_label),
            y_label=copy_method(self.y_label),
            z_label=copy_method(self.z_label),
            parameters=copy_method(self.parameters),
            process_history=copy_method(self.process_history)
        )

    @property
    def number_of_signals(self):
        return len(self.y)

    def pop(self, index: int) -> Signal:
        sig = self.get_signal(index)
        self.delete(index)
        return sig

    def delete(self, index: int | Iterable[int] | slice):
        if isinstance(index, int):
            index = [index]
        index.sort(reverse=True)  # delete largest to smallest to avoid issue of changing index
        for i in index:
            self.z = np.delete(self.z, i, axis=0)
            self.y = np.delete(self.y, i)

    def get_signal(self, y_index: int, process_history: bool = True, copy_: bool = False) -> Signal:
        """

        Parameters
        ----------
        y_index
        process_history:
            True: get x, z
            False: get x, z
        copy_:
            True: data will be a copy.
            False: data will be a view (until edited)

        Returns
        -------

        Should return a 'view' and not 'copy'. But will become a copy if edited.
        https://numpy.org/doc/stable/user/basics.copies.html

        """
        x, y = self.x, self.z[y_index, :]

        if copy_:
            x, y = x.copy(), y.copy()

        sig = self._signal(x=x, y=y, x_label=self.x_label, y_label=self.y_label,
                           name=f"slice_{self.y_label}: {self.y[y_index]}", id_=y_index)
        sig.extract_value = self.y[y_index]
        return sig

    def signal_iter(self) -> Iterator[Signal]:
        for i in range(self.number_of_signals):
            yield self.get_signal(i, process_history=True)

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
        np.savez(path, x=self.x, y=self.y, z=self.z, name=self.name,
                 x_label=self.x_label, y_label=self.y_label, z_label=self.z_label, **kwargs
                 )

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
        x_label, y_label, z_label, name = npzfile['x_label'], npzfile['y_label'], npzfile['z_label'], npzfile['name']
        return cls(x, y, z, x_label=x_label, y_label=y_label, z_label=z_label, name=name)
