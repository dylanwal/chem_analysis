def search_by_retention_time_single(
        picking_library: PickingLibrary,
        peak: PeakContinuous,
        a_tolerance: int | float = 0.1,
        number_of_matches: int = 1,
) -> PeakCompound:
    if isinstance(peak, PeakContinuous):
        raise ValueError("Invalid 'Peak' type.")

    distance = np.abs(picking_library.retention_times - peak.x)
    index = get_top_n_matches(distance, number_of_matches)
    index = index[distance[index] < a_tolerance]

    if len(index) == 1:
        return PeakCompound(peak, picking_library.compounds[index[0]], distance[index[0]])
    elif len(index) == 0:
        return PeakCompound(peak, None, None)

    compounds = [picking_library.compounds[i] for i in index]
    distances = distance[index]
    return PeakCompound(peak, compounds, distances)