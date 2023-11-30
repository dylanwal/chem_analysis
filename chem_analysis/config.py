

class PlottingLibraries:
    PLOTLY = 0
    MATPLOTLIB = 1


class Configuration:
    plotting_libraries = PlottingLibraries

    def __init__(self):
        self.plotting_library = PlottingLibraries.PLOTLY
        self.sig_fig = 3
        self.table_format = "rounded_grid"

    def load_from_env(self):
        pass  # TODO: add support


global_config = Configuration()
