
from chem_analysis.analysis.peak_picking.picking_result import ResultPeaks, ResultPeakArray
from chem_analysis.analysis.peak_picking.library_search.picking_library import PickingLibrary


def find_peaks_retention_time_library(
        peak_result: ResultPeaks | ResultPeakArray,
        library: PickingLibrary,
) -> ResultPeaks | ResultPeakArray:
    if isinstance(peak_result, ResultPeaks):
        return find_peaks_library_single(peak_result, library)
    elif isinstance(peak_result, ResultPeakArray):
        results = ResultPeakArray(signal2d=peak_result.signal2d)
        for i in range(len(results)):
            results.add_result(find_peaks_library_single(peak_result[i], library))
        return results


def find_peaks_library_single(peak_result: ResultPeaks, library: PickingLibrary) -> ResultPeaks:
    results = ResultPeaks(peak_result.signal)

    for peak in peak_result.peaks:
        if peak in library:
            results.add_peak(peak)

    return results
