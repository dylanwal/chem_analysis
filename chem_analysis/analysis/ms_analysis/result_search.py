# from typing import Iterator, Sequence
#
# import numpy as np
#
# from chem_analysis.analysis.peak import Peak
# from chem_analysis.analysis.peak_result import ResultPeaks
# # from chem_analysis.library import Compound
#
#
# class PeakCompound:
#     def __init__(self,
#                  peak: Peak,
#                  compound: Compound | Sequence[Compound] | None,
#                  metric: int | float | np.ndarray | None,
#                  ):
#         self.peak = peak
#         self.compound = compound
#         self.metric = metric
#
#     def __str__(self):
#         return str(self.compound)
#
#     def __repr__(self):
#         return self.__str__()
#
#     def __getattr__(self, name):
#         try:
#             return getattr(self.peak, name)
#         except AttributeError:
#             return getattr(self.compound, name)
#
#     def select_one(self, compound: int | Compound):
#         """ reducing multiple compounds to one """
#         if isinstance(compound, Compound):
#             if not isinstance(self.compound, Sequence):
#                 raise ValueError(f"'{type(self).__name__}.select_one' is only for reducing multiple compounds to one.")
#             compound = self.compound.index(compound)
#
#         if isinstance(compound, int):
#             if 0 < compound < len(self.compound) - 1:
#                 raise IndexError(f"Only {len(self.compound)} compounds. Index outside range.")
#             self.compound = self.compound[compound]
#             self.metric = self.metric[compound]
#             return
#
#         raise TypeError("Not valid type.")
#
#
# class ResultCompoundSearch(ResultPeaks):
#     def __init__(self, peaks: list[PeakCompound] | None = None):
#         super().__init__(peaks)
#
#     def __iter__(self) -> Iterator[PeakCompound]:
#         return iter(self.peaks)
#
#     def add_peak(self, peak: PeakCompound):
#         self.peaks.append(peak)
