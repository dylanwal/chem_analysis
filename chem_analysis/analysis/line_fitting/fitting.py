from typing import Sequence
from itertools import chain

import numpy as np
from scipy.optimize import curve_fit

from chem_analysis.analysis.line_fitting.peak_models import PeakModel


class PeaksMultiple:
    def __init__(self, peaks: Sequence[PeakModel]):
        self.peaks = peaks

    def __call__(self, x: np.ndarray, *args) -> np.ndarray:
        self.set_args(args)

        y = np.zeros_like(x, dtype=x.dtype)
        for i, peak in enumerate(self.peaks):
            y += peak(x)

        return y

    def get_args(self) -> tuple:
        return tuple(chain(*(peak.get_args() for peak in self.peaks)))

    def set_args(self, args):
        args = list(args)
        try:
            for peak in self.peaks:
                args_ = args[:peak.number_args]
                del args[:peak.number_args]
                peak.set_args(args_)
        except IndexError:
            raise IndexError("too few args")
        if len(args) != 0:
            raise IndexError("not all args used")


class ResultPeakFitting:
    def __init__(self):
        self.multipeak: PeaksMultiple | None = None
        self.covariance = None

    @property
    def peaks(self) -> Sequence:
        return self.multipeak.peaks


# class ResultPeakStats:
#     def __init__(self, parent: list[PeakStats]):
#         self.parent = parent
#
#     @property
#     def tallest_peak(self) -> None | PeakModel:
#         if len(self.peaks) == 0:
#             return None
#
#         return max(self.peaks, key=lambda x: x.max_value)
#
#     def stats(self) -> list[OrderedDict]:
#         out = []
#         for peak in self.peaks:
#             out.append(peak.stats())
#
#         return out
#
#     def print_stats(self, sig_figs: int = 3, output_str: bool = False, **kwargs):
#         """ Prints stats out for peak. """
#         from tabulate import tabulate
#
#         if "tablefmt" not in kwargs:
#             kwargs["tablefmt"] = "simple_grid"
#
#         stats_list = self.stats()
#         rows = []
#         for stats in stats_list:
#             cols = []
#             for value in stats.values():
#                 if isinstance(value, float) or isinstance(value, int):
#                     value = apply_sig_figs(value, sig_figs)
#                 cols.append(value)
#             rows.append(cols)
#
#         text = tabulate(rows, stats_list[0].keys(), **kwargs)
#         if output_str:
#             return text
#
#         print(text)


def peak_deconvolution(
        peaks: Sequence[PeakModel],
        xdata: np.ndarray,
        ydata: np.ndarray,
        **kwargs
) -> ResultPeakFitting:
    result = ResultPeakFitting()
    multipeak = PeaksMultiple(peaks)

    output = curve_fit(
        f=multipeak,
        xdata=xdata,
        ydata=ydata,
        p0=multipeak.get_args(),
        **kwargs
    )
    args, covariance = output
    multipeak.set_args(args)

    result.multipeak = multipeak
    result.covariance = covariance
    return result


def peak_deconvolution_with_n_peaks(
        peaks: Sequence[type],
        xdata: np.ndarray,
        ydata: np.ndarray,
        **kwargs
) -> ResultPeakFitting:
    peaks = [peak() for peak in peaks]

    # get initial conditions


    result = ResultPeakFitting()
    multipeak = PeaksMultiple(peaks)

    output = curve_fit(
        f=multipeak,
        xdata=xdata,
        ydata=ydata,
        p0=multipeak.get_args(),
        **kwargs
    )
    args, covariance = output
    multipeak.set_args(args)

    result.multipeak = multipeak
    result.covariance = covariance
    return result


def peak_deconvolution_auto(
        peaks: Sequence[type],
        xdata: np.ndarray,
        ydata: np.ndarray,
        **kwargs
) -> ResultPeakFitting:
    # do magic

    # get initial conditions

    result = ResultPeakFitting()
    multipeak = PeaksMultiple(peaks)

    output = curve_fit(
        f=multipeak,
        xdata=xdata,
        ydata=ydata,
        p0=multipeak.get_args(),
        **kwargs
    )
    args, covariance = output
    multipeak.set_args(args)

    result.multipeak = multipeak
    result.covariance = covariance
    return result
