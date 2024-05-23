from collections import OrderedDict

import numpy as np
import plotly.graph_objs as go

import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library.JSON"
chemistry_lib = ca.mass_spec.GCLibrary.from_JSON(lib_path)

def get_time_series(folder_path: str, pattern: str) -> ca.gc_lc.GCMSSignal2D:
    import glob

    files = glob.glob(folder_path + "/" + pattern)
    files.sort(key=lambda x: int(x[66:-9]))
    times = np.array([int(f[66:-9]) for f in files])
    data = [ca.gc_lc.GCParser.from_Agilent_D_folder(f)[0] for f in files]

    return ca.gc_lc.GCMSSignal2D.from_signals(data, times)


def first_load(folder_: str, pattern: str):
    data = get_time_series(folder_, pattern)
    print(data)
    data.to_npy(folder_ + r"\data")


def process_single(signal):
    signal.processor.add(ca.processing.baseline.SectionMinMax(sections=100, window=15, number_of_deviations=4))

    # signal_ = data.get_signal(3)
    # fig = go.Figure()
    # fig.add_trace(go.Scatter(x=signal_.x, y=signal_.y))
    # fig.show()

    peaks = ca.analysis.peak_picking.find_peaks_scipy(signal, scipy_kwargs={"height": 20000, "width": 0.1})
    peak_result = ca.analysis.boundary_detection.rolling_ball(peaks, n=10, min_height=0.05, n_points_with_pos_slope=1)

    picking_lib = chemistry_lib.to_picking_library()
    peaks_compounds = ca.analysis.peak_picking.library_search.find_peaks_retention_time_library(peak_result, picking_lib)
    # print(peaks_compounds)
    # print(peaks_compounds.stats_table().to_csv_str())

    return peaks_compounds


class ResultCompound:
    def __init__(self, parent, compound: ca.gc_lc.Compound, area: np.ndarray):
        self.parent = parent
        self.compound = compound
        self.area = area

    @property
    def mmol(self) -> float:
        return self.area/self.parent.calibrate.area*self.parent.mmol

    def add(self, i: int, area: float):
        self.area[i] = area


class ResultGrouper:
    def __init__(self, time_: np.ndarray, mmol: float = 1):
        self.time = time_
        self.compounds: list[ResultCompound] = []
        self.calibrate = None
        self.mmol = mmol

    def __getitem__(self, item: int):
        return self.compounds[item]

    @property
    def groups(self) -> list[str]:
        return list({comp.compound.group for comp in self.compounds})

    def get_by_group(self, group: str) -> list[ResultCompound]:
        return [comp for comp in self.compounds if comp.compound.group == group]

    def add(self, i: int, compound: ca.analysis.peak_picking.picking_peak.PeakCompound):
        for comp in self.compounds:
            if comp.compound is compound.compound:
                comp.add(i, compound.stats.area)
                return

        new_comp = ResultCompound(self, compound.compound, np.zeros_like(self.time))
        new_comp.add(i, compound.stats.area)
        self.compounds.append(new_comp)

    def add_result(self, i: int, result: chem_analysis.analysis.picking_result.ResultPeaks):
        for peak in result:
            self.add(i, peak)

    def set_calibrate(self, compound: ca.gc_lc.Compound):
        for comp in self.compounds:
            if comp.compound is compound:
                self.calibrate = comp
                return
        raise ValueError("No calibrat found")


def map_with_fill_zeros(times: list[float], comp_times: list[float], values: list[float]) -> list[float]:
    out = []
    for t in times:
        if t in comp_times:
            index = comp_times.index(t)
            out.append(values[index])
        else:
            out.append(0)
    return out


def lighter_shades(color: str, n: int) -> list[str]:
    hex_symbol_flag = False
    if color[0] == "#":
        hex_symbol_flag = True
        color = color[1:]

    # Convert hex color to RGB
    r = int(color[0:2], 16)
    g = int(color[2:4], 16)
    b = int(color[4:6], 16)

    # Calculate step size for each color channel
    r_step = (255 - r) / n
    g_step = (255 - g) / n
    b_step = (255 - b) / n

    # Generate lighter shades
    shades = []
    for i in range(n):
        r_new = int(r + i * r_step)
        g_new = int(g + i * g_step)
        b_new = int(b + i * b_step)
        new_color = "{:02x}{:02x}{:02x}".format(r_new, g_new, b_new)
        if hex_symbol_flag:
            new_color = "#" + new_color
        shades.append(new_color)

    return shades


class ResultPlotter:
    # picked with HSB (H:variable, Saturation: 100, Brightness:30)
    COLORS = OrderedDict()
    COLORS["substrate"] = "000000"
    COLORS["misc"] = "545454"
    COLORS["peroxide"] = "545400"
    COLORS["alcohol TMS"] = "540000"
    COLORS["alcohol acetylated"] = "633500"
    COLORS["ketone"] = "046104"
    COLORS["carboxylic acid TMS"] = "000054"
    COLORS["hydroxy acids TMS"] = "540054"
    COLORS["hydroxy acids acetylated TMS"] = "540054"
    COLORS["dicarboxylic acid TMS"] = "003030"

    @classmethod
    def add_products(cls, fig: go.Figure, results):
        groups = results.groups

        lengend_names = []
        for group_name, group_color in reversed(cls.COLORS.items()):
            if group_name not in groups:
                continue
            group = groups[group_name]
            group.sort(key=lambda x: x.compound.smiles.molar_mass)
            colors = lighter_shades(group_color, len(group))

            for compound, color in zip(group, colors):
                if compound.group not in lengend_names:
                    lengend_names.append(compound.group)
                    kwargs = dict(name=compound.group)
                else:
                    kwargs = dict(showlegend=False)

                kwargs.update(dict(
                        hoverinfo="text",
                        text=compound.compound.label,
                        mode="lines",
                        stackgroup="one",
                        line=dict(color=f'#{color}'),
                        legendgroup=compound.compound.group,
                        hoveron='points+fills',
                ))

                mmols = map_with_fill_zeros(results.times, compound.times, compound.mmols)

                fig.add_trace(go.Scatter(x=results.times, y=mmols, **kwargs))

    @classmethod
    def add_substrate(cls, fig: go.Figure, results, substrate_mmol: int | float = None):
        substrate = results.substrate
        fig.add_trace(go.Scatter(x=substrate.times, y=substrate.mmols, mode="lines", name=substrate.compound.name,
                                 line=dict(color="black"), legendgroup=substrate.compound.name))

        if substrate_mmol is not None:
            fig.add_trace(go.Scatter(
                x=[0, results.times[-1]],
                y=[substrate_mmol, substrate_mmol],
                mode="lines", line=dict(color="gray", dash="dash"), showlegend=False, legendgroup=substrate.compound.name))

            fig.add_trace(go.Scatter(x=results.times, y=[result.total_mmol for result in results.results],
                                 mode="lines", line=dict(color="gray"), name="mass balance", legendgroup="decane"))

    @classmethod
    def add_layout(cls,
                   fig: go.Figure,
                   x_axis: str = "",
                   y_axis: str = "",
                   title: str = ""
                   ):
        fig.layout.font = dict(family="Arial", size=18, color="black")
        fig.layout.height = 600
        fig.layout.width = 1200
        fig.layout.plot_bgcolor = "white"

        fig.layout.xaxis.title = x_axis
        fig.layout.xaxis.tickprefix = "<b>"
        fig.layout.xaxis.ticksuffix = "</b>"
        fig.layout.xaxis.showline = True
        fig.layout.xaxis.linewidth = 5
        fig.layout.xaxis.mirror = True
        fig.layout.xaxis.linecolor = "black"
        fig.layout.xaxis.showgrid = False
        # fig.layout.xaxis.range = [0, 1000]

        fig.layout.yaxis.title = y_axis
        fig.layout.yaxis.tickprefix = "<b>"
        fig.layout.yaxis.ticksuffix = "</b>"
        fig.layout.yaxis.showline = True
        fig.layout.yaxis.linewidth = 5
        fig.layout.yaxis.mirror = True
        fig.layout.yaxis.linecolor = "black"
        fig.layout.yaxis.showgrid = False
        # fig.layout.yaxis.range = [0, 12]

        fig.layout.title.text = title


def main(folder_: str):
    data = ca.gc_lc.GCMSSignal2D.from_file(folder_ + r"\data.npz")
    results = ResultGrouper(time_=data.y)
    for sig in range(len(data)):
        results.add_result(sig, process_single(data.get_signal(sig)))
    TCB = chemistry_lib.find_by_label("TCB")
    results.set_calibrate(TCB)

    fig = go.Figure(layout={"template": ca.global_config.plotly_layout()})
    total = np.zeros_like(results.calibrate.mmol)
    for comp in results.compounds:
        if comp.compound is TCB:
            continue
        fig.add_trace(go.Scatter(x=results.time, y=comp.mmol, name=comp.compound.label))
        total += comp.mmol
    fig.add_trace(go.Scatter(x=results.time, y=total, name="mass balance", line={"dash": "dash"}))
    fig.add_trace(go.Scatter(x=[results.time[0], results.time[-1]], y=[total[0], total[0]], line={"color": "black"}, showlegend=False))
    fig.layout.xaxis.title = "<b>time (min)</b>"
    fig.layout.yaxis.title = "<b>mmol</b>"
    fig.layout.yaxis.range = [0, None]
    fig.show()


    fig = go.Figure()
    groups = results.groups
    lengend_names = []
    for group_name, group_color in reversed(ResultPlotter.COLORS.items()):
        if group_name not in groups:
            continue
        group = results.get_by_group(group_name)
        group.sort(key=lambda x: len(str(x.compound.smiles)))
        colors = lighter_shades(group_color, len(group))

        for compound, color in zip(group, colors):
            if compound.compound.group not in lengend_names:
                lengend_names.append(compound.compound.group)
                kwargs = dict(name=compound.compound.group)
            else:
                kwargs = dict(showlegend=False)

            kwargs.update(dict(
                    hoverinfo="text",
                    text=compound.compound.label,
                    mode="lines",
                    stackgroup="one",
                    line=dict(color=f'#{color}'),
                    legendgroup=compound.compound.group,
                    hoveron='points+fills',
            ))

            fig.add_trace(go.Scatter(x=results.time, y=compound.mmol, **kwargs))

    ResultPlotter.add_layout(fig,
                         x_axis="<b>Time (min) </b>",
                         y_axis="<b>mmol of chemical</b>",
                         title="<b>MC of Decane (100 psi O<sub>2</sub> flow-through, 160 <sup>o</sup>C) </b>",
                         )
    fig.show()


if __name__ == "__main__":
    folder = r"C:\Users\nicep\Desktop\research_wis\data\11\11_23\gc_ms"
    # first_load(folder, "DJW-11-23-*min-TMS.D")
    main(folder)
