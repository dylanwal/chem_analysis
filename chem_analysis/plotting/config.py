import enum


class NormalizationOptions(enum.Enum):
    NONE = 0
    AREA = 1
    PEAK_HEIGHT = 2


class SignalColorOptions(enum.Enum):
    SINGLE = 0
    DIVERSE = 1


class PlotConfig:
    NORMALIZATION_OPTIONS = NormalizationOptions
    COLOR_OPTIONS = SignalColorOptions

    def __init__(self):
        # signal
        self.normalize: NormalizationOptions = NormalizationOptions.NONE

        # peak
        # self.peak_show_trace: bool = False
        self.peak_show_shade: bool = True
        self.peak_show_bounds: bool = False
        self.peak_show_max: bool = True
        self.signal_connect_gaps: bool = False
        self.signal_color: str | SignalColorOptions | None = None
        self.peak_bound_color: str | None = None
        self.peak_bound_line_width: float = 3
        self.peak_bound_height: float = 0.05  # % of max
        self.peak_marker_size: float = 3
