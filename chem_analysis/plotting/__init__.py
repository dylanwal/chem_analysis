from chem_analysis import global_config
from chem_analysis.plotting.plotting import signal, add_peaks
from chem_analysis.plotting.plotting2D import (signal2D_contour, signal2D_surface,
                                               signal2D_overlap_signals, signal2D_slices_peaks, signal2D_stack_signals)

# if global_config.PLOTTING_LIBRARIES.PLOTLY in global_config.get_plotting_options():
#     import chem_analysis.plotting.plotly_plots.plotly_utils as plotly_utils
# if global_config.PLOTTING_LIBRARIES.MATPLOTLIB in global_config.get_plotting_options():
#     from plotly_plots.plotly_config import PlotlyConfig
# if global_config.PLOTTING_LIBRARIES.PYGRAPHQT in global_config.get_plotting_options():
#     from plotly_plots.plotly_config import PlotlyConfig
