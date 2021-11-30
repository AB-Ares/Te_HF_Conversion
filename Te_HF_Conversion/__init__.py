"""
Te_HF_Conversion
============================
Te_HF_Conversion provides a function and an example script
for converting a given elastic thickness of the lithosphere
into a planetary heat flow, or a temperature profile
into a yield strength envelope and elastic thickness. 
Radiogenic heat contribution and wet or dry rheologies are 
implemented.

   Conversion_Te_HF
      Determine the surface, crust, and mantle heat flows
      given the input parameters.

   Curv_Moment
      Performs the  to moment calculation.

   Conversion_Tprofile_Te
      Convert a temperature profile into an
      elastic thickness and yield strength envelope.

   Brittle_Strength
      Compute the brittle part of the yield strength
      envelope based on planetary constants.

   Planet_constants
      Contains the constants used in the calculations.
"""
from ._version import get_versions

from .CalcYSE_HF import Conversion_Te_HF
from .CalcYSE_HF import Conversion_Tprofile_Te
from .CalcYSE_HF import Curv_Moment
from .CalcYSE_HF import Brittle_Strength
from .Planet_constants import *

del CalcYSE_HF

__version__ = get_versions()["version"]
del get_versions

__author__ = "Adrien Broquet"

__all__ = [
    "Conversion_Te_HF",
    "Curv_Moment",
    "Conversion_Tprofile_Te",
    "Brittle_Strength",
    "Planet_constants",
]
