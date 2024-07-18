import enum


class NormalizationOptions(enum.Enum):
    NONE = 0
    AREA = 1
    PEAK_HEIGHT = 2


class SignalColorOptions(enum.Enum):
    SINGLE = 0
    DIVERSE = 1
