import pathlib
from typing import Sequence

import numpy as np

import chem_analysis.utils.math as utils_math
from chem_analysis.library.library import Library
from chem_analysis.library.chemical import Chemical
from chem_analysis.library.attribute import MassSpectrum
from chem_analysis.library.condition import Method


class MSLibrary:
    def __init__(self,
                 chemicals: Sequence[Chemical],
                 ms_mass: np.ndarray | None = None,
                 ms_intensity: np.ndarray | None = None,
                 ):
        """
        >> compounds and retention times need to be reordered such that mass_spec is a continuous array. This is done for
        efficiency reasons in algorithms. if a compound has a mass_spec and no retention time_set retention time to -1.

        Parameters
        ----------
        chemicals: Iterable[Chemical]
            compounds
        ms_mass:
            m/z values
        ms_intensity:
            can be size less than compounds
            normalize to total ion count

        """
        self.chemicals = chemicals
        self.ms_mass = ms_mass
        self.ms_intensity = ms_intensity

    def __iter__(self):
        return iter(self.chemicals)

    def __getitem__(self, item: int | slice) -> Chemical:
        return self.chemicals[item]

    def to_npz(self, path: str | pathlib.Path, **kwargs):
        chemicals = [chem.name for chem in self.chemicals]
        np.savez(path, ms_mass=self.ms_mass, ms_intensity=self.ms_intensity, chemicals=chemicals, **kwargs)

    @classmethod
    def from_npz(cls, path: str | pathlib.Path, **kwargs):
        npzfile = np.load(str(path), allow_pickle=True, **kwargs)
        ms_mass, ms_intensity, chemicals = npzfile['ms_mass'], npzfile['ms_intensity'], npzfile['chemicals']
        return cls(chemicals, ms_mass=ms_mass, ms_intensity=ms_intensity)

    @classmethod
    def from_library(cls,
                     library: Library,
                     method: str | None = None,
                     grab_any_ms: bool | str | tuple[str] = False,
                     ms_dtype: str | np.dtype = np.uint16,
                     ):
        chemicals = []
        mass_specs = []
        for chem in library:
            ms = None
            for attr in chem.attributes:
                if isinstance(attr, MassSpectrum):
                    if not grab_any_ms and method is not None:
                        if any(cond.value == method for cond in attr.conditions if isinstance(cond, Method)):
                            ms = attr.value
                    else:
                        if ms is None:
                            ms = attr.value

            if ms is not None:
                mass_specs.append(ms)
                chemicals.append(chem)

        if len(chemicals) == 0:
            raise ValueError("No mass_specs found and thus no Library was built.")

        # sort for ms to make algorithms efficient
        none_count = sum(True for ms in mass_specs if ms is None)
        if none_count < len(mass_specs):
            if none_count >= 1:
                mass_specs, chemicals, retention_times = move_nones_to_end(mass_specs, chemicals)
                first_none = mass_specs.index(None)
                ms_mass, ms_intensity = unify_ms_values(mass_specs[:first_none], ms_dtype)
            else:
                ms_mass, ms_intensity = unify_ms_values(mass_specs, ms_dtype)
        else:
            ms_mass, ms_intensity = None, None

        return cls(chemicals, ms_mass, ms_intensity)


def move_nones_to_end(sort_list, *other_lists):
    filtered_tuples = [t for t in zip(sort_list, *other_lists) if t[0] is not None]
    none_tuples = [t for t in zip(sort_list, *other_lists) if t[0] is None]
    combined_tuples = filtered_tuples + none_tuples
    separated_lists = list(zip(*combined_tuples))

    return [list(lst) for lst in separated_lists]


def unify_ms_values(mass_specs: list[np.ndarray], dtype: np.dtype = np.uint16) -> tuple[np.ndarray, np.ndarray]:
    # x are not the same
    max_x = np.max([np.max(ms[:, 0]) for ms in mass_specs])
    min_x = np.min([np.min(ms[:, 0]) for ms in mass_specs])

    x = np.arange(min_x, max_x+1, dtype=utils_math.min_int_dtype(max_x, min_x))
    y = np.zeros((len(mass_specs), len(x)), dtype=dtype)
    for i, sig in enumerate(mass_specs):
        y_new = utils_math.map_discrete_x_axis(x, np.round(sig[:, 0]), sig[:, 1])
        y[i, :] = (y_new/np.max(y_new) * np.iinfo(dtype).max).astype(dtype) # normalize to total ion count

    return x, y
