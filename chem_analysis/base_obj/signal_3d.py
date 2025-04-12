from typing import Sequence, Iterable
import pathlib
import copy

import numpy as np

from chem_analysis.base_obj.parameters import Parameters
import chem_analysis.utils.math as math_utils
from chem_analysis.base_obj.unify_methods_2d import UnifyMethod2D, UnifyMethodStrict2D
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.base_obj.signal_ import check_array_inf_nan


def validate_input(x: np.ndarray, y: np.ndarray, z: np.ndarray, w: np.ndarray):
    if len(x.shape) != 1:
        raise ValueError(f"'x' must shape 1. \n\treceived: {x.shape}")
    if len(y.shape) != 1:
        raise ValueError(f"'y' must shape 1. \n\treceived: {y.shape}")
    if len(z.shape) != 1:
        raise ValueError(f"'z' must shape 1. \n\treceived: {z.shape}")
    if len(w.shape) != 3:
        raise ValueError(f"'w' must shape 3. \n\treceived: {w.shape}")
    if x.shape[0] != w.shape[2]:
        raise ValueError(f"'x' and 'w[2]' must have same shape. \n\treceived: x:{x.shape} "
                         f"|| w.shape[2]: {w.shape[2]}")
    if y.shape[0] != w.shape[1]:
        raise ValueError(f"'y' and 'w[1]' must have same shape. \n\treceived: y:{y.shape} "
                         f"|| w.shape[1]: {w.shape[1]}")
    if z.shape[0] != w.shape[0]:
        raise ValueError(f"'z' and 'w[0]' must have same shape. \n\treceived: z:{z.shape} "
                         f"|| w.shape[0]: {w.shape[0]}")
    check_array_inf_nan(x, name='x')
    check_array_inf_nan(y, name='y')
    check_array_inf_nan(z, name='z')
    check_array_inf_nan(w, name='w')


class Signal3D:
    """ signal 3D

    A signal is any x-y-z-w data.

    """
    _signal = Signal2D
    __count = 0

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
        parameters: Parameters
            various meta-data
        process_history: list[str]
            List of processing methods the signal has been through
        """
        validate_input(x, y, z, w)

        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.id_ = id_ or Signal3D.__count
        Signal3D.__count += 1
        self.name = name or f"signal3D_{self.id_}"
        self.x_label = x_label or "x_axis"
        self.y_label = y_label or "y_axis"
        self.z_label = z_label or "z_axis"
        self.w_label = w_label or "w_axis"

        self.parameters = parameters
        self.extract_value = None

        if isinstance(process_history, str):
            process_history = [process_history]
        self.process_history = process_history or []

    def __repr__(self):
        text = f"{self.name}: "
        text += f"{self.x_label} vs {self.y_label} vs {self.z_label} vs {self.w_label}"
        text += f" (shape: {self.w.shape})"
        return text

    def copy_with(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, w: np.ndarray, deep: bool = True):
        if deep:
            copy_method = copy.deepcopy
        else:
            copy_method = copy.copy

        return Signal3D(
            x,
            y,
            z,
            w,
            name=copy_method(self.name),
            x_label=copy_method(self.x_label),
            y_label=copy_method(self.y_label),
            z_label=copy_method(self.z_label),
            w_label=copy_method(self.w_label),
            parameters=copy_method(self.parameters),
            process_history=copy_method(self.process_history)
        )

    def pop(self, index: int) -> Signal2D:
        sig = self.get_signal(index)
        self.delete(index)
        return sig

    def delete(self, index: int | Iterable):
        if isinstance(index, int):
            index = [index]
        index.sort(reverse=True)  # delete largest to smallest to avoid issue of changing index
        for i in index:
            self.w = np.delete(self.w, i, axis=0)
            self.z = np.delete(self.z, i)

    def get_signal(self, z_index: int, process_history: bool = True, copy_: bool = False) -> Signal2D:
        """

        Parameters
        ----------
        z_index
        process_history:
            True: get x, y, w
            False: get x, y, w
        copy_:
            True: data will be a copy.
            False: data will be a view (until edited)

        Returns
        -------

        Should return a 'view' and not 'copy'. But will become a copy if edited.
        https://numpy.org/doc/stable/user/basics.copies.html

        """
        x, y, z = self.x, self.y, self.w[z_index, :]
        if copy_:
            x, y, z = x.copy(), y.copy(), z.copy()

        sig = self._signal(x=x, y=y, z=z, x_label=self.x_label, y_label=self.y_label, z_label=self.z_label,
                           name=f"slice_{self.z_label}: {self.z[z_index]}", id_=z_index)
        sig.extract_value = self.z[z_index]
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
    def write_npz(self, path: str | pathlib.Path, sparse: bool = False, **kwargs):
        """Save an array to a binary file in NumPy ``.npz`` format."""
        if sparse:
            from chem_analysis.utils.sparse_data import numpy_to_sparse
            coords, data, shape = numpy_to_sparse(self.w)
            coords.astype(math_utils.min_uint_dtype(np.max(coords)))
            np.savez(path, x=self.x, y=self.y, z=self.z, coords=coords,
                     x_label=self.x_label, y_label=self.y_label, z_label=self.z_label, w_label=self.w_label,
                     data=data, shape=shape, **kwargs)
        else:
            np.savez(path, x=self.x, y=self.y, z=self.z, w=self.w,
                     x_label=self.x_label, y_label=self.y_label, z_label=self.z_label, w_label=self.w_label,
                     **kwargs)

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

        npzfile = np.load(str(path))
        x, y, z = npzfile['x'], npzfile['y'], npzfile['z']
        x_label, y_label, z_label, w_label, name = npzfile['x_label'], npzfile['y_label'], npzfile['z_label'], npzfile[
            'w_label'], npzfile['name']
        return cls(x, y, z, w, x_label=x_label, y_label=y_label, z_label=z_label, w_label=w_label, name=name)
