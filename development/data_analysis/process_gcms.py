import pathlib

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library_color.json"
LIBRARY = ca.library.Library.from_JSON(lib_path)
picking_lib_fid = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_fid")
picking_lib_ms = ca.analysis.ms_analysis.PickingLibrary.from_library(LIBRARY, "decane_ms")

parameters = ""


class CompoundTimeSeries:
    def __init__(self, compound: ca.library.Compound):
        self.compound = compound
        self.areas = []
        self.times = []

    def __getattr__(self, name):
        return getattr(self.compound, name)

    def add_time_point(self, time_: int | float, area: int | float):
        self.times.append(time_)
        self.areas.append(area)


class ResultTimeSeries:
    def __init__(self):
        self.compounds: list[CompoundTimeSeries] = []
        self._times: list[int | float] = []
        self.compounds_by_times: list[list[CompoundTimeSeries]] = []
        self._compound_list: list[ca.gc_lc.Compound] = []

    @property
    def times(self) -> np.ndarray:
        return np.array(self._times)

    def add_result(self, result: ca.analysis.ms_analysis.ResultCompoundSearch, time_: int | float):
        comp_for_time_ = []
        for comp in result.peaks:
            if comp.compound is not None:
                comp_out = self._add_compound(comp, time_)
                comp_for_time_.append(comp_out)

        self.compounds_by_times.append(comp_for_time_)
        self._times.append(time_)

    def _add_compound(self,
                      compound: ca.analysis.ms_analysis.PeakCompound,
                      time_: int | float,
                      ) -> CompoundTimeSeries:
        if compound.compound in self._compound_list:
            index = self._compound_list.index(compound.compound)
            comp = self.compounds[index]
            comp.add_time_point(time_, compound.peak.area())
        else:
            self._compound_list.append(compound.compound)
            comp = CompoundTimeSeries(compound.compound)
            comp.add_time_point(time_, compound.peak.area())
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
            compounds.append(comp.compound.label)
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


def process_single(data_path: pathlib.Path, label: str, zero):
    # load data
    ms_file = data_path / f"{label}_data.ms"
    fid_file = data_path / f"{label}_FID1A.ch"
    ini_file = data_path / f"{label}_pre_post.ini"
    ms, fid = ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)

    # ms_file = data_path / f"{zero}_data.ms"
    # fid_file = data_path / f"{zero}_FID1A.ch"
    # ini_file = data_path / f"{zero}_pre_post.ini"
    # ms_zero, fid_zero = ca.gc_lc.GCParser.from_Agilent_D_files(ini_file, ms_file, fid_file)
    # span_ = ca.p.edit.ReplaceSpans(value=0, x_spans=((None, 3.7), (5.5, None)))
    # x, y = span_.run(fid_zero.x_raw, fid_zero.y_raw)

    # fid
    fid.processor.add(
        ca.p.edit.ReplaceSpans(value=0,
                               x_spans=(
                                   (1.3, 3.7),  # solvent
                                   # (44.15, 44.85)  # PPh3
                               )
                               ),
        # ca.p.baseline.SubtractOptimize(x=x, y=y),
        ca.p.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4)
    )
    peak_locations = ca.a.peak_picking.find_peaks_scipy(fid, scipy_kwargs={"height": 500, "width": 0.15})
    peaks = ca.a.integration.rolling_ball(peak_locations, n=5, min_height=0.002, n_points_with_pos_slope=2)
    fid_compounds = ca.a.ms_analysis.search_by_retention_time(picking_lib_fid, peaks)

    fid_fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    ca.plotting.signal(fid, fig=fid_fig)
    ca.plotting.peaks(fid_compounds, fig=fid_fig)
    fid_fig.layout.title = "FID"

    # ms
    ms.processor.add(
        # ca.p.edit.ReplaceSpans(value=0,
        #                        x_spans=(
        #                            (44.1, 45.1)  # PPh3
        #                        )
        #                        ),
        ca.p.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4)
    )
    peak_locations = ca.a.peak_picking.find_peaks_scipy(ms, scipy_kwargs={"height": 8000, "width": 0.1})
    peaks = ca.a.integration.rolling_ball(peak_locations, n=5, min_height=0.002, n_points_with_pos_slope=2)
    ms_compounds = ca.a.ms_analysis.search_by_retention_time(picking_lib_ms, peaks)

    ms_fig = go.Figure(layout=ca.plotting.plotly_utils.layout())
    ca.plotting.signal(ms, fig=ms_fig)
    ca.plotting.peaks(ms_compounds, fig=ms_fig)
    ms_fig.layout.title = "MS"

    print("finished analyzing:", label)
    global parameters
    parameters += "\n" + str(data_path) + "\n\t ms:" + str(ms.processor.methods) + "\n\t fid:" + str(
        ms.processor.methods)
    return ms_fig, fid_fig, ms_compounds, fid_compounds


def process_timeseries(data_path: pathlib.Path, data_label: str, labels: list[str], times: np.ndarray):
    global parameters
    parameters += str(data_path) + "\n" + str(data_label) + "\n"
    # process data
    fid_compounds, ms_compounds, fid_figs, ms_figs = [], [], [], []
    for label in labels:
        ms_fig_, fid_fig_, ms_comp, fid_comp = process_single(data_path, label, zero=labels[0])
        ms_figs.append(ms_fig_)
        fid_figs.append(fid_fig_)
        ms_compounds.append(ms_comp)
        fid_compounds.append(fid_comp)

    # re-organizing data
    fid_timeseries = ResultTimeSeries()
    ms_timeseries = ResultTimeSeries()
    for i in range(len(times)):
        fid_timeseries.add_result(fid_compounds[i], times[i])
        ms_timeseries.add_result(ms_compounds[i], times[i])

    # saving data
    ca.plotting.plotly_utils.merge_figures(fid_figs, filename=data_path / (data_label + '_fid.html'))
    ca.plotting.plotly_utils.merge_figures(ms_figs, filename=data_path / (data_label + '_ms.html'))
    ms_timeseries.to_csv(data_path / (data_label + '_ms.csv'))
    data = fid_timeseries.to_csv(data_path / (data_label + '_fid.csv'))
    print(data)
    with open(data_path / (data_label + "_params.txt"), mode='w') as f:
        f.write(parameters)


def main():
    data_path = pathlib.Path(rf"C:\Users\nicep\Desktop\research_wis\data\11\11_52\GCMS")
    times = np.array([0, 10, 20, 30, 45, 60, 120])
    data_label = "DJW-11-52-V2"
    labels = [data_label + f"-t{i}" for i in times]
    process_timeseries(data_path, data_label, labels, times)


if __name__ == "__main__":
    main()
