from chem_analysis import global_config
from chem_analysis.plotting.plotting import signal, baseline, calibration, peaks

if global_config.PLOTTING_LIBRARIES.PLOTLY in global_config.get_plotting_options():
    from chem_analysis.plotting.plotly_plots.plotly_config import PlotlyConfig
# if global_config.PLOTTING_LIBRARIES.MATPLOTLIB in global_config.get_plotting_options():
#     from plotly_plots.plotly_config import PlotlyConfig
# if global_config.PLOTTING_LIBRARIES.PYGRAPHQT in global_config.get_plotting_options():
#     from plotly_plots.plotly_config import PlotlyConfig
