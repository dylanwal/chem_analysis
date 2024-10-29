
import numpy as np

from chem_analysis.mass_spec import MSSignal
import chem_analysis.utils.math as utils_math


def unify_ms_signals(signal1: MSSignal, signal2: MSSignal) -> tuple[MSSignal, MSSignal]:
    if not np.all(np.isclose(signal1.x, signal2.x, rtol=0.01)):
        # unify mz
        max_mz = int(max(signal1.mz.max(), signal2.mz.max()))
        x = np.arange(max_mz, dtype=utils_math.min_int_dtype(max_mz, 0))
        y_1 = utils_math.map_discrete_x_axis(x, np.round(signal1.mz), signal1.y)
        y_2 = utils_math.map_discrete_x_axis(x, np.round(signal2.mz), signal2.y)

    return signal1, signal2
