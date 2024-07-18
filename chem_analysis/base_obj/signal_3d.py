from typing import Sequence, Iterable
import pathlib

import numpy as np

import chem_analysis.utils.math as math_utils
from chem_analysis.base_obj.unify_methods_2d import UnifyMethod2D, UnifyMethodStrict2D
from chem_analysis.processing.processor import Processor
from chem_analysis.analysis.peak import PeakParent2D
from chem_analysis.base_obj.signal_2d import Signal2D


def validate_input(x_raw: np.ndarray, y_raw: np.ndarray, z_raw: np.ndarray, w_raw: np.ndarray):
    if len(x_raw.shape) != 1:
        raise ValueError(f"'x_raw' must shape 1. \n\treceived: {x_raw.shape}")
    if len(y_raw.shape) != 1:
        raise ValueError(f"'y_raw' must shape 1. \n\treceived: {y_raw.shape}")
    if len(z_raw.shape) != 1:
        raise ValueError(f"'z_raw' must shape 1. \n\treceived: {z_raw.shape}")
    if len(w_raw.shape) != 3:
        raise ValueError(f"'w_raw' must shape 3. \n\treceived: {w_raw.shape}")
    if x_raw.shape[0] != w_raw.shape[2]:
        raise ValueError(f"'x_raw' and 'w_raw[2]' must have same shape. \n\treceived: x_raw:{x_raw.shape} "
                         f"|| w_raw.shape[2]: {w_raw.shape[2]}")
    if y_raw.shape[0] != w_raw.shape[1]:
        raise ValueError(f"'y_raw' and 'w_raw[1]' must have same shape. \n\treceived: y_raw:{y_raw.shape} "
                         f"|| w_raw.shape[1]: {w_raw.shape[1]}")
    if z_raw.shape[0] != w_raw.shape[0]:
        raise ValueError(f"'z_raw' and 'w_raw[0]' must have same shape. \n\treceived: z_raw:{z_raw.shape} "
                         f"|| w_raw.shape[0]: {w_raw.shape[0]}")


class Signal3D:
    """ signal 3D

    A signal is any x-y-z-w data.

    """
    _signal = Signal2D
    __count = 0
    _peak_type = PeakParent2D  # TODO: upgrade to 2D bounded

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 w: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 w_label: str = None,
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
        z: np.ndarray[k]
            raw z data, length k
        w: np.ndarray[k,j,i]
            raw z data, shape k,j,i
        x_label: str
            x-axis label
        y_label: str
            y-axis label
        z_label: str
            z-axis label
        w_label: str
            4-axis label
        name: str
            user defined name
        """
        validate_input(x, y, z, w)

        self.x_raw = x
        self.y_raw = y
        self.z_raw = z
        self.w_raw = w
        self.id_ = id_ or Signal3D.__count
        Signal3D.__count += 1
        self.name = name or f"signal3D_{self.id_}"
        self.x_label = x_label or "x_axis"
        self.y_label = y_label or "y_axis"
        self.z_label = z_label or "z_axis"
        self.w_label = w_label or "w_axis"

        self.processor = Processor()
        self._x = None
        self._y = None
        self._z = None
        self._w = None

        self.extract_value = None

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label} vs {self.z_label} vs {self.w_label}"
        text += f" (shape: {self.w_raw.shape})"
        return text

    def _process(self):
        self._x, self._y, self._z, self._w = self.processor.run(self.x_raw, self.y_raw, self.z_raw, self.w_raw)

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
    def w(self) -> np.ndarray:
        if not self.processor.processed:
            self._process()
        return self._w

    def pop(self, index: int) -> Signal2D:
        sig = self.get_signal(index)
        self.delete(index)
        return sig

    def delete(self, index: int | Iterable):
        if isinstance(index, int):
            index = [index]
        index.sort(reverse=True)  # delete largest to smallest to avoid issue of changing index
        for i in index:
            self.w_raw = np.delete(self.w_raw, i, axis=0)
            self.z_raw = np.delete(self.z_raw, i)

    def get_signal(self, z_index: int, processed: bool = True, copy_: bool = False) -> Signal2D:
        """

        Parameters
        ----------
        z_index
        processed:
            True: get x, y, w
            False: get x_raw, y_raw, w_raw
        copy_:
            True: data will be a copy.
            False: data will be a view (until edited)

        Returns
        -------

        Should return a 'view' and not 'copy'. But will become a copy if edited.
        https://numpy.org/doc/stable/user/basics.copies.html

        """
        if processed:
            x, y, z = self.x, self.y, self.w[z_index, :]
        else:
            x, y, z = self.x_raw, self.y_raw, self.w_raw[z_index, :]

        if copy_:
            x, y, z = x.copy(), y.copy(), z.copy()

        sig = self._signal(x=x, y=y, z=z, x_label=self.x_label, y_label=self.y_label,  z_label=self.z_label,
                           name=f"slice_{self.z_label}: {self.z[z_index]}", id_=z_index)
        sig.extract_value = self.z[z_index]
        if not processed:
            sig.processor = self.processor.get_copy()

        return sig

    @classmethod
    def from_signals(cls,
                     signals: Sequence[Signal2D],
                     z: np.ndarray = None,
                     z_label: str = None,
                     unify_method: UnifyMethod2D = UnifyMethodStrict2D(),
                     ):  # -> Signal3D
        """ Turn Sequence of Signal2Ds into a Signal3D"""
        if z is None:
            z = np.arange(len(signals))
        else:
            if len(z.shape) != 1 and z.shape[0] == len(signals):
                raise ValueError("The number of signals must be the same as the number of z points.\n"
                                 f"\tnumber of signals: {len(signals)}\n\tnumber of z points:{z.shape[0]}")

        x_label = signals[0].x_label
        y_label = signals[0].y_label
        w_label = signals[0].z_label

        x, y, w = unify_method.run(signals)
        return cls(x, y, z, w, x_label=x_label, y_label=y_label, z_label=z_label, w_label=w_label)

    ####################################################################################################################
    ## Save/Load from file #############################################################################################
    ####################################################################################################################
    def to_npz(self, path: str | pathlib.Path, sparse: bool = False, **kwargs):
        """Save an array to a binary file in NumPy ``.npz`` format."""
        if sparse:
            from chem_analysis.utils.sparse_data import numpy_to_sparse
            coords, data, shape = numpy_to_sparse(self.w)
            coords.astype(math_utils.min_uint_dtype(np.max(coords)))
            np.savez(path, x=self.x, y=self.y, z=self.z, coords=coords, data=data, shape=shape, **kwargs)
        else:
            np.savez(path, x=self.x, y=self.y, z=self.z, w=self.w, **kwargs)

    @classmethod
    def from_npz(cls, path: str | pathlib.Path):
        npzfile = np.load(str(path))
        if 'coords' in npzfile.files:
            # sparse array
            from chem_analysis.utils.sparse_data import sparse_to_numpy
            coords = npzfile['coords']
            data = npzfile['data']
            shape = npzfile['shape']
            w = sparse_to_numpy(coords, data, shape)
        else:
            w = npzfile['w']

        x, y, z = npzfile['x'], npzfile['y'], npzfile['z']
        x_label = y_label = z_label = w_label = None
        return cls(x, y, z, w, x_label=x_label, y_label=y_label, z_label=z_label, w_label=w_label)

