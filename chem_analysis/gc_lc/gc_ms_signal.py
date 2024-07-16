import copy

import numpy as np

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.gc_lc.gc_parameters import GCParameters
from chem_analysis.mass_spec.ms_signal_2D import MSSignal2D


class GCMSSignal(Signal):
    """
    Gas Chromatogram Signal
    """
    def __init__(self,
                 ms_raw: MSSignal2D,
                 x_label: str = None,
                 y_label: str = None,
                 parameters: GCParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "retention time"
        y_label = y_label or "intensity"
        name = name or ms_raw.name
        y_raw = np.sum(ms_raw.z_raw, axis=1)
        super().__init__(copy.copy(ms_raw.y), y_raw, x_label, y_label, name, id_)
        self.parameters = parameters
        self.ms_raw = ms_raw
