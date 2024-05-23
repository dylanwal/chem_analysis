

def to_picking_library(library, method: str):
    from chem_analysis.analysis.peak_picking.library_search.picking_library import PickingLibrary
    from chem_analysis.analysis.peak_picking.picking_peak import PeakForPickingCompound
    peaks = [PeakForPickingCompound(compound) for compound in library if method in compound.methods]
    return PickingLibrary(peaks)

# def get_n_nearest_compounds(self, retention_time: float, n: int = 1) -> list[Compound]:
#     distance = []
#     for i, peak in enumerate(self.peaks):
#         distance[i] = abs(peak.retention_time - retention_time)
#
#     sort_index = np.argsort(distance)
#     compounds = []
#     for i in range(n):
#         compounds.append(self.compounds[sort_index[i]])
#     return compounds