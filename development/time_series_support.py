from collections import OrderedDict

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca


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
            raise ValueError("Internal standard not found in result.")

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


def map_with_fill_zeros(times: list[float], comp_times: list[float], values: list[float]) -> list[float]:
    out = []
    for t in times:
        if t in comp_times:
            index = comp_times.index(t)
            out.append(values[index])
        else:
            out.append(0)
    return out


def plot_results(results: ResultTimeSeries, groups: list[str], *, fig: go.Figure = None) -> go.Figure:
    if fig is None:
        fig = go.Figure()

    compounds = results.compounds

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

            fig.add_trace(go.Scatter(x=compound.times, y=compound.mmols, **kwargs))

    return fig

# def add_substrate(cls, fig: go.Figure, results: ResultTimeSeries, substrate_mmol: int | float = None):
#     substrate = results.substrate
#     fig.add_trace(go.Scatter(x=substrate.times, y=substrate.mmols, mode="lines", name=substrate.compound.name,
#                              line=dict(color="black"), legendgroup=substrate.compound.name))
#
#     if substrate_mmol is not None:
#         fig.add_trace(go.Scatter(
#             x=[0, results.times[-1]],
#             y=[substrate_mmol, substrate_mmol],
#             mode="lines", line=dict(color="gray", dash="dash"), showlegend=False, legendgroup=substrate.compound.name))
#
#         fig.add_trace(go.Scatter(x=results.times, y=[result.total_mmol for result in results.results],
#                              mode="lines", line=dict(color="gray"), name="mass balance", legendgroup="decane"))
#
