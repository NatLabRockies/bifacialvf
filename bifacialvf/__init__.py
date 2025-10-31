# ensure python3 compatible division and printing
from __future__ import division, print_function, absolute_import
from bifacialvf.bifacialvf import simulate, skycomposition_method, getEPW, readInputTMY  # main program
from bifacialvf.vf import getBackSurfaceIrradiances, getFrontSurfaceIrradiances, getGroundShadeFactors  # main subroutines
from bifacialvf.vf import getSkyConfigurationFactors, trackingBFvaluescalculator, rowSpacing # helper functions
from bifacialvf.sun import hrSolarPos, perezComp, solarPos, sunIncident # solar position and value
from bifacialvf.loadVFresults import loadVFresults # utility for reading result files
from bifacialvf.BF_BifacialIrradiances.LandscapeSingleHour import LandscapeSingleHour # For calculateBilInterpol
from bifacialvf.BF_BifacialIrradiances.PortraitSingleHour import PortraitSingleHour # For calculateBilInterpol

try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:
    # for python < 3.8 (remove when dropping 3.7 support)
    from importlib_metadata import PackageNotFoundError, version

try:
    __version__ = version(__package__)
except PackageNotFoundError:
    __version__ = "0+unknown"