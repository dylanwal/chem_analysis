import logging
from typing import Sequence

import numpy as np

from chem_analysis.base_obj.unify_methods_2d import UnifyMethod2D, UnifyMethodMS2D
from chem_analysis.base_obj.signal_2d import Signal2D
from chem_analysis.mass_spec.ms_signal_3D import MSSignal3D
from chem_analysis.gc_lc.gc_ms_signal import GCMSSignal

logger = logging.getLogger(__name__)


class GCMSSignal2D(Signal2D):
    _signal = GCMSSignal

    def __init__(self,
                 ms: MSSignal3D,
                 x_label: str = None,
                 y_label: str = None,
                 z_label: str = None,
                 name: str = None
                 ):
        x_label = x_label or "retention time"
        y_label = y_label or "time"
        z_label = z_label or "intensity"
        z_raw = np.sum(ms.w_raw, axis=2)
        super().__init__(ms.y_raw, ms.z_raw, z_raw, x_label, y_label, z_label, name)
        self.ms = ms

    def get_signal(self, y_index: int, processed: bool = False) -> GCMSSignal:
        return super().get_signal(y_index, processed)

    @classmethod
    def from_signals(cls,
                     signals: Sequence[GCMSSignal],
                     y: np.ndarray = None,
                     y_label: str = None,
                     unify_method: UnifyMethod2D = UnifyMethodMS2D(),
                     ):  # -> Signal2D
        if y is None:
            y = np.arange(len(signals))
        else:
            if len(y.shape) != 1 and y.shape[0] == len(signals):
                raise ValueError("The number of signals must be the same as the number of y points.\n"
                                 f"\tnumber of signals: {len(signals)}\n\tnumber of y points:{y.shape[0]}")

        x_label = signals[0].x_label
        z_label = signals[0].y_label
        ms = MSSignal3D.from_signals([sig.ms_raw for sig in signals], y, y_label, unify_method)
        return cls(ms=ms, x_label=x_label, y_label=y_label, z_label=z_label)

