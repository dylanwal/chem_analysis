from typing import Sequence

import numpy as np

import chem_analysis.utils.math as utils_math
from chem_analysis.library import Library, Compound


class PickingLibrary:
    def __init__(self,
                 compounds: Sequence[Compound],
                 retention_times: np.ndarray | None = None,
                 ms_mass: np.ndarray | None = None,
                 ms_intensity: np.ndarray | None = None,
                 ):
        """
        >> compounds and retention times need to be reordered such that mass_spec is a continuous array. This is done for
        efficiency reasons in algorithms. if a compound has a mass_spec and no retention time_set retention time to -1.

        Parameters
        ----------
        compounds: Iterable[Compound]
            compounds
        retention_times:
            can be size less than compounds
            -1 means no retention time
        ms_mass:
            m/z values
        ms_intensity:
            can be size less than compounds
            normalize to total ion count

        """
        self.compounds = compounds
        self.retention_times = retention_times
        self.ms_mass = ms_mass
        self.ms_intensity = ms_intensity

    def __iter__(self):
        return iter(self.compounds)

    def __getitem__(self, item: int | slice) -> Compound:
        return self.compounds[item]

    @classmethod
    def from_library(cls,
                     library: Library,
                     method: str,
                     grab_any_ms: bool | str | tuple[str] = False,
                     ms_dtype: str | np.dtype = np.uint16,
                     ):
        compounds = []
        retention_times = []
        mass_specs = []
        for compound in library:
            if method in compound.methods:
                compounds.append(compound)
                response = compound.get_response(method)
                mass_specs.append(response.mass_spectrum)
                retention_times.append(response.retention_time)
            if grab_any_ms and mass_specs[-1] is None:
                ms = compound.get_ms()
                if ms is not None:
                    mass_specs[-1] = ms

        count_None = sum([1 for i in mass_specs if i is None])
        if count_None == len(mass_specs):
            mass_specs = None
        if all(retention_times) is None:
            retention_times = None
        if mass_specs is None and retention_times is None:
            raise ValueError("No mass_specs or retention_times found and thus no PickingLibrary was built.")

        # sort for ms to make algorithms efficient
        if count_None > 0 and mass_specs is not None:
            if None in mass_specs:
                mass_specs, compounds, retention_times = move_nones_to_end(mass_specs, compounds, retention_times)
                first_none = mass_specs.index(None)
                ms_mass, ms_intensity = unify_ms_values(mass_specs[:first_none], ms_dtype)
            else:
                ms_mass, ms_intensity = unify_ms_values(mass_specs, ms_dtype)
        else:
            ms_mass, ms_intensity = None, None

        return PickingLibrary(compounds, np.array(retention_times), ms_mass, ms_intensity)


def move_nones_to_end(sort_list, *other_lists):
    filtered_tuples = [t for t in zip(sort_list, *other_lists) if t[0] is not None]
    none_tuples = [t for t in zip(sort_list, *other_lists) if t[0] is None]
    combined_tuples = filtered_tuples + none_tuples
    separated_lists = list(zip(*combined_tuples))

    return [list(lst) for lst in separated_lists]


def unify_ms_values(mass_specs: list[np.ndarray], dtype: np.dtype = np.uint16) -> tuple[np.ndarray, np.ndarray]:
    # x are not the same
    max_x = np.max([np.max(ms[0]) for ms in mass_specs])
    min_x = np.min([np.min(ms[0]) for ms in mass_specs])

    x = np.arange(min_x, max_x+1, dtype=utils_math.min_int_dtype(max_x, min_x))
    y = np.zeros((len(mass_specs), len(x)), dtype=dtype)
    for i, sig in enumerate(mass_specs):
        y_new = utils_math.map_discrete_x_axis(x, np.round(sig[0]), sig[1])
        y[i, :] = (y_new/np.max(y_new) * np.iinfo(dtype).max).astype(dtype) # normalize to total ion count

    return x, y
