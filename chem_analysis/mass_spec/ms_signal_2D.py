from __future__ import annotations
from typing import Sequence

import numpy as np

from chem_analysis.base_obj.unify_methods import UnifyMethod, UnifyMethodMS
from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.mass_spec.ms_signal import MSSignal
from chem_analysis.mass_spec.ms_parameters import MSParameters


class MSSignal2D(Signal2D):
    _signal = MSSignal

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 parameters: MSParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "time"
        z_label = z_label or "counts"
        super().__init__(x, y, z, x_label, y_label, z_label, name, id_)
        self.parameters = parameters

    def get_signal(self, index: int, processed: bool = True) -> MSSignal:
        return super().get_signal(index, processed)

    @classmethod
    def from_signals(cls,
                     signals: Sequence[Signal],
                     y: np.ndarray = None,
                     y_label: str = None,
                     unify_method: UnifyMethod = UnifyMethodMS()
                     ) -> MSSignal2D:
        return super().from_signals(signals, y, y_label, unify_method)
