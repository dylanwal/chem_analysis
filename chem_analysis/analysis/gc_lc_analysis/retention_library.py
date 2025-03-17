import pathlib
from typing import Sequence

import numpy as np

from chem_analysis.library.library import Library
from chem_analysis.library.chemical import Chemical
from chem_analysis.library.attribute import RetentionTime
from chem_analysis.library.condition import Method


class RetentionTimeLibrary:
    def __init__(self,
                 chemicals: Sequence[Chemical],
                 times: np.ndarray | None = None,
                 ):
        """
        retention time

        Parameters
        ----------
        chemicals: Iterable[Chemical]
            compounds
        times:
            m/z values

        """
        self.chemicals = chemicals
        self.times = times

    def __iter__(self):
        return iter(self.chemicals)

    def __getitem__(self, item: int | slice) -> Chemical:
        return self.chemicals[item]

    def to_npz(self, path: str | pathlib.Path, **kwargs):
        chemicals = [chem.name for chem in self.chemicals]
        np.savez(path, times=self.times, chemicals=chemicals, **kwargs)

    @classmethod
    def from_npz(cls, path: str | pathlib.Path, **kwargs):
        npzfile = np.load(str(path), allow_pickle=True, **kwargs)
        times, chemicals = npzfile['times'], npzfile['chemicals']
        return cls(chemicals, times=times)

    @classmethod
    def from_library(cls,
                     library: Library,
                     method: str | None = None,
                     ):
        chemicals = []
        times = []
        for chem in library:
            time_ = None
            for attr in chem.attributes:
                if isinstance(attr, RetentionTime):
                    if method is None:
                        time_ = attr.value
                    elif any(cond.value == method for cond in attr.conditions if isinstance(cond, Method)):
                        time_ = attr.value

            if time_ is not None:
                times.append(time_)
                chemicals.append(chem)

        if len(chemicals) == 0:
            raise ValueError("No retention_times found and thus no Library was built.")

        # sort
        times = np.array(times)
        index = np.argsort(times)
        times = times[index]
        chemicals_sorted = [chemicals[i] for i in index]

        return cls(chemicals_sorted, times)
