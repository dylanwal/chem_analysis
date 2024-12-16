import pathlib
from typing import Sequence

import numpy as np

import chem_analysis as ca
from chem_analysis.analysis.peak_result import ResultPeaks
from chem_analysis.analysis.peak import PeakContinuous


class TimeSeries:
    def __init__(self, times: Sequence[float], compounds: Sequence[str], values: np.ndarray):
        self.times = times
        self.compounds = compounds
        self.values = values

    def to_csv(self, filename: pathlib.Path = None) -> str:
        text = ""
        times = np.insert(self.times, 0, 0)
        text = ",".join([str(i) for i in times]) + "\n"
        values = self.values.T
        for i in range(len(self.compounds)):
            text += ",".join([self.compounds[i]] + [str(i_) for i_ in values[i, :]]) + "\n"

        if filename is not None:
            # TODO: improve file handling/naming
            with open(filename, "w") as f:
                f.write(text)
        return text


def results_to_timeseries(data: Sequence[ResultPeaks], times: Sequence[float]) -> TimeSeries:
    if len(times) != len(data):
        raise ValueError(f"Length of times does not match length of data: "
                         f"\n\tlen(data):{len(data)}"
                         f"\n\tlen(times):{len(times)}"
                         )

    result = ResultTimeSeries()
    for result_, t in zip(data, times):
        result.add_result(result_, t)

    compounds, times, areas = result.to_numpy()
    return TimeSeries(times, compounds, areas)


class CompoundTimeSeries:
    def __init__(self, compound: ca.library.Compound):
        self.compound = compound
        self.areas = []
        self.times = []

    def __getattr__(self, name: str):
        if hasattr(self.compound, name):
            return getattr(self.compound, name)
        return self.compound

    def add_time_point(self, time_: int | float, area: int | float):
        self.times.append(time_)
        self.areas.append(area)


class ResultTimeSeries:
    def __init__(self):
        self.compounds: list[CompoundTimeSeries] = []
        self._times: list[int | float] = []
        self.compounds_by_times: list[list[CompoundTimeSeries]] = []
        self._compound_list: list = []

    @property
    def times(self) -> np.ndarray:
        return np.array(self._times)

    def add_result(self, result: ResultPeaks, time_: int | float):
        comp_for_time_ = []
        for comp in result.peaks:
            if comp.label is not None:
                comp_out = self._add_compound(comp, time_)
                comp_for_time_.append(comp_out)

        self.compounds_by_times.append(comp_for_time_)
        self._times.append(time_)

    def _add_compound(self,
                      compound: PeakContinuous,
                      time_: int | float,
                      ) -> CompoundTimeSeries:
        if compound.label in self._compound_list:
            index = self._compound_list.index(compound.label)
            comp = self.compounds[index]
            comp.add_time_point(time_, compound.properties.area())
        else:
            self._compound_list.append(compound.label)
            comp = CompoundTimeSeries(compound.label)
            comp.add_time_point(time_, compound.properties.area())
            self.compounds.append(comp)

        return comp

    def to_numpy(self) -> tuple[list[str], np.ndarray, np.ndarray]:
        times = set()
        for comp in self.compounds:
            for t in comp.times:
                times.add(t)

        times = np.array(list(times))
        times.sort()

        compounds = []
        areas = np.zeros((len(times), len(self.compounds)))
        for i, comp in enumerate(self.compounds):
            compounds.append(str(comp.label))
            for ii, t in enumerate(comp.times):
                index = np.argmin(abs(times - t))
                areas[index, i] = comp.areas[ii]

        return compounds, times, areas

    def to_csv(self, filename: pathlib.Path = None) -> str:
        compounds, times, mmols = self.to_numpy()

        text = ""
        times = np.insert(times, 0, 0)
        text = ",".join([str(i) for i in times]) + "\n"
        mmols = mmols.T
        for i in range(len(compounds)):
            text += ",".join([compounds[i]] + [str(i_) for i_ in mmols[i, :]]) + "\n"

        if filename is not None:
            # TODO: improve file handling/naming
            with open(filename, "w") as f:
                f.write(text)
        return text
