from src.chem_analysis.analysis.algorithms.despike import despike
from src.chem_analysis.analysis.algorithms.baseline import poly_baseline


despike_methods = {"default": despike}
baseline_methods = {"polynomial": poly_baseline}
smoothing_methods = {}