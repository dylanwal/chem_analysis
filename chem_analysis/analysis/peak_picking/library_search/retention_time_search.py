
from chem_analysis.analysis.peak_picking.picking_result import ResultPeakPicking, ResultPeakPickingArray
from chem_analysis.analysis.peak_picking.library_search.picking_library import PickingLibrary


def find_peaks_retention_time_library(
        peak_result: ResultPeakPicking | ResultPeakPickingArray,
        library: PickingLibrary,
) -> ResultPeakPicking | ResultPeakPickingArray:
    if isinstance(peak_result, ResultPeakPicking):
        return find_peaks_library_single(peak_result, library)
    elif isinstance(peak_result, ResultPeakPickingArray):
        results = ResultPeakPickingArray(signal2d=peak_result.signal2d)
        for i in range(len(results)):
            results.add_result(find_peaks_library_single(peak_result[i], library))
        return results


def find_peaks_library_single(peak_result: ResultPeakPicking, library: PickingLibrary) -> ResultPeakPicking:
    results = ResultPeakPicking(peak_result.signal)

    for peak in peak_result.peaks:
        if peak in library:
            results.add_peak(peak)

    return results
