import logging
import copy
from collections import OrderedDict
from pathlib import Path

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca

logger = logging.getLogger(__name__)


class CompoundTimeSeries:
    def __init__(self, compound: ca.gc_lc.Compound, method: str):
        self.method = method
        self.compound = compound
        self._areas = []
        self._times = []
        self._internal_areas = []
        self._internal_mmol = []

    def __getattr__(self, name):
        return getattr(self.compound, name)

    def add_time_point(self,
                       time_: int | float,
                       area: int | float,
                       internal_area: int | float,
                       internal_mmol: int | float
                       ):
        self._times.append(time_)
        self._areas.append(area)
        self._internal_areas.append(internal_area)
        self._internal_mmol.append(internal_mmol)

    @property
    def times(self) -> np.ndarray:
        return np.array(self._times)

    @property
    def areas(self) -> np.ndarray:
        return np.array(self._areas)

    @property
    def areas_normalized(self) -> np.ndarray:
        return self.areas/np.array(self._internal_areas)

    @property
    def mmols(self) -> np.ndarray:
        response_factor = self.compound.get_response(self.method).response
        return self.areas_normalized * response_factor * self._internal_mmol

    @property
    def groups(self) -> list[str]:
        return self.compound.groups


class ResultTimeSeries:
    def __init__(self,
                 method: str,
                 internal_standard: ca.gc_lc.Compound,
                 internal_standard_mmol: int | float = 1
                 ):
        self.method = method
        self.internal_standard = internal_standard
        self.internal_standard_mmol = internal_standard_mmol
        self.compounds: list[CompoundTimeSeries] = []
        self._times: list[int | float] = []
        self.compounds_by_times: list[list[CompoundTimeSeries]] = []
        self._compound_list: list[ca.gc_lc.Compound] = []

    @property
    def times(self) -> np.ndarray:
        return np.array(self._times)

    def add_result(self, result: ca.analysis.ms_analysis.ResultCompoundSearch, time_: int | float):
        for comp in result.peaks:
            if comp.compound is self.internal_standard:
                internal_standard = comp
                break
        else:
            logger.error(f"Internal standard not found in result ({result.name}). Results dropped")
            return

        comp_for_time_ = []
        for comp in result.peaks:
            if comp is internal_standard:
                continue
            comp_out = self._add_compound(comp, time_, internal_standard.peak.area())
            comp_for_time_.append(comp_out)

        self.compounds_by_times.append(comp_for_time_)
        self._times.append(time_)

    def _add_compound(self,
                      compound: ca.analysis.ms_analysis.PeakCompound,
                      time_: int | float,
                      internal_standard_area: int | float
                      ) -> CompoundTimeSeries:
        if compound.compound is None:
            return  #TODO:

        if compound.compound in self._compound_list:
            index = self._compound_list.index(compound.compound)
            comp = self.compounds[index]
            comp.add_time_point(time_, compound.peak.area(), internal_standard_area, self.internal_standard_mmol)
        else:
            self._compound_list.append(compound.compound)
            comp = CompoundTimeSeries(compound.compound, self.method)
            comp.add_time_point(time_, compound.peak.area(), internal_standard_area, self.internal_standard_mmol)
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
        mmols = np.zeros((len(times), len(self.compounds)))
        for i, comp in enumerate(self.compounds):
            compounds.append(comp.compound.label)
            for ii, t in enumerate(comp.times):
                index = np.argmin(abs(times-t))
                mmols[index, i] = comp.mmols[ii]

        return compounds, times, mmols

    def to_csv(self, filename: str | Path = None) -> str:
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


def map_with_fill_zeros(times: list[float], comp_times: list[float], values: list[float]) -> list[float]:
    out = []
    for t in times:
        if t in comp_times:
            index = comp_times.index(t)
            out.append(values[index])
        else:
            out.append(0)
    return out


def plot_results(
        results: ResultTimeSeries,
        groups: list[str],
        *,
        fig: go.Figure = None,
        add_zero: bool = False,
        carbon: bool = False
) -> go.Figure:
    if fig is None:
        fig = go.Figure()

    compounds = copy.copy(results.compounds)

    # sort
    groups.append("misc")
    groups_sorted = OrderedDict()
    for group_name in groups:
        if len(compounds) == 0:
            break

        group_ = []
        FLAG = False
        for i in range(len(compounds)):
            for i in range(len(compounds)):
                if group_name in compounds[i].groups:
                    group_.append(compounds.pop(i))
                    break
            else:
                FLAG = True
            if FLAG:
                break

        if len(group_) == 0:
            continue
        # sort within group
        group_.sort(key=lambda x: x.compound.smiles.molar_mass)
        group_ = list(reversed(group_))
        groups_sorted[group_name] = group_

    lengend_names = []
    for group_name, group_ in groups_sorted.items():
        for compound in group_:
            if group_name not in lengend_names:
                lengend_names.append(group_name)
                kwargs = dict(name=group_name)
            else:
                kwargs = dict(showlegend=False)

            kwargs.update(dict(
                    hoverinfo="text",
                    text=compound.compound.label,
                    mode="lines",
                    stackgroup="one",
                    line=dict(color=compound.color),
                    fillcolor=compound.color,
                    legendgroup=group_name,
                    hoveron='points+fills',
            ))

            # mmols = map_with_fill_zeros(results.times, compound.times, compound.mmols)


            if add_zero:
                if compound.times[0] == 0:
                    x = compound.times
                    y = compound.mmols
                elif compound.compound.label == "C10":
                    x = np.insert(compound.times, 0, 0)
                    y = np.insert(compound.mmols, 0, 0.0658)
                else:
                    x = np.insert(compound.times, 0, 0)
                    y = np.insert(compound.mmols, 0, 0)
            else:
                x = compound.times
                y = compound.mmols

            if carbon:
                y = y * (str(compound.compound.smiles).count("C") - str(compound.compound.smiles).count("[Si]")*3)

            fig.add_trace(go.Scatter(x=x, y=y, **kwargs))

    return fig
