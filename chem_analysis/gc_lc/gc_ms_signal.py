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
                 x: np.ndarray = None,
                 y: np.ndarray = None,
                 ms: MSSignal2D = None,
                 x_label: str = None,
                 y_label: str = None,
                 parameters: GCParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        if x is None or y is None:
            if ms is None:
                raise ValueError('x and y must not be None or ms must not be None')
            x = copy.copy(ms.y)
            y = np.sum(ms.z, axis=1)

        x_label = x_label or "retention time"
        y_label = y_label or "intensity"
        name = name or ms.name
        super().__init__(x, y, x_label, y_label, name, id_)
        self.parameters = parameters
        self.ms_raw = ms
