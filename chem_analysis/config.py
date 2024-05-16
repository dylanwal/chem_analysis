import enum


class PlottingLibraries(enum.Enum):
    PLOTLY = 0
    MATPLOTLIB = 1
    PYGRAPHQT = 2


class Configuration:
    plotting_libraries = PlottingLibraries

    def __init__(self):
        self.preferred_plot = PlottingLibraries.PLOTLY
        self._plotting_libraries = []
        self._find_available_plotting_libraries()

        self.sig_fig: int = 3
        self.table_format: str = "rounded_grid"
        self.processing_save_intermediates: bool = False
        self.max_mz: int = 1000

    def load_from_env(self):
        pass  # TODO: add support

    def get_plotting_options(self) -> list[PlottingLibraries]:
        self._find_available_plotting_libraries()
        if self._plotting_libraries is None:
            raise RuntimeError("No plotting libraries installed. Please install one of the following:"
                               "\n\tplotly: `pip install plotly'"
                               "\n\tmatplotlib: 'pip install matplotlib'"
                               "\n\tpygraphqt: 'pip install pygraphqt'")

        if self.preferred_plot in self._plotting_libraries:
            self._plotting_libraries.remove(self.preferred_plot)
            self._plotting_libraries.insert(0, self.preferred_plot)

        return self._plotting_libraries

    def _find_available_plotting_libraries(self):
        try:
            import plotly
            self._plotting_libraries.append(PlottingLibraries.PLOTLY)
        except ImportError:
            pass
        try:
            import matplotlib
            self._plotting_libraries.append(PlottingLibraries.MATPLOTLIB)
        except ImportError:
            pass
        try:
            import pyqtgraph
            self._plotting_libraries.append(PlottingLibraries.PYGRAPHQT)
        except ImportError:
            pass

    @staticmethod
    def plotly_layout():
        import plotly.graph_objs as go
        template = go.layout.Template()
        template.layout.font = dict(family="Arial", size=18, color="black")
        template.layout.plot_bgcolor = "white"
        template.layout.width, template.layout.height = 1200, 600
        template.layout.xaxis.tickprefix = "<b>"
        template.layout.xaxis.ticksuffix = "<b>"
        template.layout.xaxis.showline = True
        template.layout.xaxis.linewidth = 5
        template.layout.xaxis.linecolor = "black"
        template.layout.xaxis.ticks = "outside"
        template.layout.xaxis.tickwidth = 4
        template.layout.xaxis.showgrid = False
        template.layout.xaxis.mirror = True
        template.layout.yaxis.tickprefix = "<b>"
        template.layout.yaxis.ticksuffix = "<b>"
        template.layout.yaxis.showline = True
        template.layout.yaxis.linewidth = 5
        template.layout.yaxis.linecolor = "black"
        template.layout.yaxis.ticks = "outside"
        template.layout.yaxis.tickwidth = 4
        template.layout.yaxis.showgrid = False
        template.layout.yaxis.mirror = True

        return template


global_config = Configuration()
