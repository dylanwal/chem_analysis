from chem_analysis import global_config
from chem_analysis.plotting.plotting import signal, baseline, calibration, peaks
from chem_analysis.plotting.plotting2D import (signal2D_contour, signal2D_surface,
                                               signal2D_slices, signal2D_slices_peaks, signal2D_slices_separate)

if global_config.PLOTTING_LIBRARIES.PLOTLY in global_config.get_plotting_options():
    import chem_analysis.plotting.plotly_plots.plotly_utils as plotly_utils
# if global_config.PLOTTING_LIBRARIES.MATPLOTLIB in global_config.get_plotting_options():
#     from plotly_plots.plotly_config import PlotlyConfig
# if global_config.PLOTTING_LIBRARIES.PYGRAPHQT in global_config.get_plotting_options():
#     from plotly_plots.plotly_config import PlotlyConfig
