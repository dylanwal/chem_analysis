import numpy as np

from chem_analysis.base_obj.signal_array import SignalArray


class IRArray(SignalArray):
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 z: np.ndarray,
                 x_label: str = "wave_number",
                 y_label: str = "time",
                 z_label: str = "Absorbance",
                 name: str = None
                 ):
        super().__init__(x, y, z, x_label, y_label, z_label, name)

    @property
    def names(self):
        return [i.name for i in self.signals]

    @property
    def y_labels(self):
        return [i.y_label for i in self.signals]

    @property
    def x_label(self):
        return self.signals[0].x_label

    @property
    def number_of_signals(self):
        return len(self.signals)
