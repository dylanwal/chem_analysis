import enum


class PlottingLibrary(enum.Enum):
    plotly = 0
    matplotlib = 1
    pygraphqt = 2


class Configuration:
    PLOTTING_LIBRARIES = PlottingLibrary

    def __init__(self):
        self._preferred_plot = PlottingLibrary.plotly
        self._plotting_libraries: list[PlottingLibrary] = []
        self._find_available_plotting_libraries()

        self.sig_fig: int = 4
        self.table_format: str = "rounded_grid"
        self.processing_save_intermediates: bool = False
        self.max_mz: int = 1000

    @property
    def preferred_plot(self):
        return self._preferred_plot

    @preferred_plot.setter
    def preferred_plot(self, value: str | PlottingLibrary):
        if isinstance(value, str):
            value = [i for i in PlottingLibrary if i.name == value.lower()]
            if len(value) == 0:
                raise ValueError(f"Plotting library not recognized: {value}.from "
                                 f"Try one of the following: {[f'{i.name}' for i in PlottingLibrary]}.")
            value = value[0]

        self._preferred_plot = value

    def load_from_env(self):
        pass  # TODO: add support    this should include plot config too

    def get_plotting_options(self) -> list[PlottingLibrary]:
        self._find_available_plotting_libraries()
        if self._plotting_libraries is None:
            raise RuntimeError("No plotting libraries installed. Please install one of the following:"
                               "\n\tplotly: `pip install plotly'"
                               "\n\tmatplotlib: 'pip install matplotlib'"
                               "\n\tpygraphqt: 'pip install pygraphqt'")
        return self._plotting_libraries

    def get_plotting_lib(self) -> PlottingLibrary:
        self._find_available_plotting_libraries()
        if self._plotting_libraries is None:
            raise RuntimeError("No plotting libraries installed. Please install one of the following:"
                               "\n\tplotly: `pip install plotly'"
                               "\n\tmatplotlib: 'pip install matplotlib'"
                               "\n\tpygraphqt: 'pip install pygraphqt'")

        if self._preferred_plot in self._plotting_libraries:
            return self.preferred_plot

        return self._plotting_libraries[0]

    def get_plotting_lib_from_fig(self, fig) -> PlottingLibrary:
        for i in self._plotting_libraries:
            if i is PlottingLibrary.plotly:
                import plotly.graph_objects as go
                if isinstance(fig, go.Figure):
                    return i
            if i is PlottingLibrary.matplotlib:
                import matplotlib.pyplot as plt
                if isinstance(fig, plt.Figure):
                    return i
            if i is PlottingLibrary.pygraphqt:
                import pyqtgraph
                if isinstance(fig, pyqtgraph.PlotWidget): #GraphicsLayoutWidget

                    return i
        raise ValueError(f"No plotting library found for figure: {type(fig)}")

    def _find_available_plotting_libraries(self):
        if not self._plotting_libraries:
            for plot_method in PlottingLibrary:
                try:
                    __import__(plot_method.name)
                    self._plotting_libraries.append(plot_method)
                except ImportError:
                    pass


global_config = Configuration()
