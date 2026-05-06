import copy
from typing import Sequence

import numpy as np

from chem_analysis.mass_spec import MSSignal


def unify_mz(args: Sequence[MSSignal], fill_to_zero: bool = True) -> list[MSSignal]:
    """

    Parameters
    ----------
    args: MSSignal
    fill_to_zero: bool
        True: force mz to start at zero
        False: mz starts at min(mz) which may or may not be zero

    Returns
    -------
    mz: np.ndarray [m]
    intensity: np.ndarray [m, len(args)]
    """
    mz, intensities = unify_mz_xy(*(s.to_numpy() for s in args), fill_to_zero=fill_to_zero, group=True)

    signals = []
    for i, arg in enumerate(args):
        sig = copy.copy(arg)
        sig.x = mz
        sig.y = intensities[i, :]
        signals.append(sig)

    return signals


def unify_mz_xy(*args: np.ndarray, fill_to_zero: bool = True, group: bool = False) \
        -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray]:
    """

    Parameters
    ----------
    args: np.ndarray [n, 2]
        column 1: mz values
        column 2: intensities
        accepts as many ars as you want
    fill_to_zero: bool
        True: force mz to start at zero
        False: mz starts at min(mz) which may or may not be zero
    group:
        True: intensity will be return as a single numpy array
        False: intensity will be returned as individual numpy array

    Returns
    -------
    mz: np.ndarray [m]
    intensity: np.ndarray [len(args), m]
    """
    min_mz = 0 if fill_to_zero else int(np.min([np.min(x[:, 0]) for x in args]))
    max_mz = int(np.max([np.max(x[:, 0]) for x in args]))
    mz = np.arange(min_mz, max_mz+1, dtype=args[0].dtype)

    intensity = np.zeros((len(args), len(mz)), dtype=float)
    for i, arg in enumerate(args):
        index = arg[:, 0].astype("uint")
        if min_mz != 0:
            index = index - min_mz
        intensity[i, index] = arg[:, 1]

    if group:
        return mz, intensity

    return tuple(np.column_stack((mz, row)) for row in intensity)
