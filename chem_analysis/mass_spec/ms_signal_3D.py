from __future__ import annotations
from typing import Sequence

import numpy as np

from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.base_obj.signal_3d import Signal3D
from chem_analysis.mass_spec.ms_signal_2D import MSSignal2D
from chem_analysis.mass_spec.ms_parameters import MSParameters


class MSSignal3D(Signal3D):
    _signal = MSSignal2D

    def __init__(self,
                 x_raw: np.ndarray,
                 y_raw: np.ndarray,
                 t_raw: np.ndarray,
                 data_raw: np.ndarray,
                 x_label: str = None,
                 y_label: str = None,
                 t_label: str = None,
                 z_label: str = None,
                 parameters: MSParameters = None,
                 name: str = None,
                 id_: int = None
                 ):
        x_label = x_label or "mass-to-charge"
        y_label = y_label or "retention time"
        z_label = z_label or "time"
        t_label = t_label or "counts"
        super().__init__(x_raw, y_raw, t_raw, data_raw, x_label, y_label, t_label, z_label, name, id_)
        self.parameters = parameters

    def get_signal(self, index: int, processed: bool = False) -> MSSignal2D:
        return super().get_signal(index, processed)

    @classmethod
    def from_signals(cls,
                     signals: Sequence[Signal2D],
                     z: np.ndarray = None,
                     x_label: str = None,
                     y_label: str = None,
                     z_label: str = None,
                     w_label: str = None,
                     min_x: int = None,
                     max_x: int = None
                     ) -> MSSignal3D:
        return super().from_signals(signals, z, x_label, y_label, z_label, min_x, max_x)

    # @classmethod
    # def from_list(cls,
    #               data: Sequence[Sequence[Sequence[np.ndarray]]],
    #               y: np.ndarray = None,
    #               x_label: str = None,
    #               y_label: str = None,
    #               z_label: str = None,
    #               min_x: int = None,
    #               max_x: int = None,
    #               ) -> MSSignal2D:
    #     return super().from_list(data, y, x_label, y_label, z_label, min_x, max_x)
