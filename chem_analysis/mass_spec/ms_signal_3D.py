from __future__ import annotations
from typing import Sequence

import numpy as np

from chem_analysis.base_obj.unify_methods_2d import UnifyMethod2D, UnifyMethodMS2D
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.base_obj.signal_3d import Signal3D
from chem_analysis.mass_spec.ms_signal_2D import MSSignal2D
from chem_analysis.mass_spec.ms_parameters import MSParameters


class MSSignal3D(Signal3D):
    _signal = MSSignal2D

    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 w: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 w_label: str = None,
                 parameters: MSParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "retention time"
        z_label = z_label or "time"
        w_label = w_label or "counts"
        super().__init__(x, y, z, w, x_label, y_label, z_label, w_label, name, id_)
        self.parameters = parameters

    def get_signal(self, z_index: int, processed: bool = True, copy_: bool = False) -> MSSignal2D:
        return super().get_signal(z_index, processed, copy_)

    @classmethod
    def from_signals(cls,
                     signals: Sequence[Signal2D],
                     z: np.ndarray = None,
                     z_label: str = None,
                     unify_method: UnifyMethod2D = UnifyMethodMS2D(),
                     ) -> MSSignal3D:
        return super().from_signals(signals, z, z_label, unify_method)
