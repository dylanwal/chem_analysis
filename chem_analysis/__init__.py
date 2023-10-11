import logging

logger = logging.getLogger("chem_analysis")  # pathlib.Path(sys.argv[0]).stem

import chem_analysis.algorithms.processing as processing
import chem_analysis.algorithms.analysis as analysis

import chem_analysis.sec as sec
import chem_analysis.nmr as nmr
import chem_analysis.ir as ir
import chem_analysis.mass_spec as mass_spec
import chem_analysis.uv_vis as uv_vis
