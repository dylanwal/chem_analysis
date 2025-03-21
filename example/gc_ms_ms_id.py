


ms_lib = ca.a.ms_analysis.MSLibrary.from_library(lib)
labels = []
for m in ms_peak_data[:]:
    matches = ca.a.ms_analysis.search_by_ms(
        ms_lib,
        m,
        scorer_ms=ca.a.ms_analysis.ScorerMultiple([ca.a.ms_analysis.ScorerDot(), offset_scorer]),
        filter_ms=[ca.a.ms_analysis.FilterMinScore(0.6), ca.a.ms_analysis.FilterTopNMatches(n=1)]
    )
    if len(matches) > 0:
        labels.append(matches[0])
    else:
        labels.append(None)