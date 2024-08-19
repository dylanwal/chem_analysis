from chem_analysis.config import global_config
import logging

logger = logging.getLogger("chem_analysis")
logger.addHandler(logging.StreamHandler())

import chem_analysis.processing as processing
p = proc = processing
import chem_analysis.analysis as analysis
a = analysis

import chem_analysis.base_obj as base
import chem_analysis.sec as sec
import chem_analysis.nmr as nmr
import chem_analysis.ir as ir
import chem_analysis.gc_lc as gc_lc
import chem_analysis.mass_spec as ms
import chem_analysis.uv_vis as uv_vis
import chem_analysis.utils
import chem_analysis.plotting as plot

