"""
Te_HF_Conversion
============================
Te_HF_Conversion provides a function and an example script
for converting a given elastic thickness of the lithosphere
into a planetary heat flow. Radiogenic heat contribution and
wet or dry rheologies are implemented.

   Conversion_Te_HF
      Determine the surface, crust, and mantle heat flows
      given the input parameters.

   Curv_Moment
      Performs the  to moment calculation.

   Planet_constants
      Contains the constants used in the calculations.
"""
from ._version import get_versions

from .CalcYSE_HF import Conversion_Te_HF
from .CalcYSE_HF import Conversion_Tprofile_Te
from .CalcYSE_HF import Curv_Moment
from .Planet_constants import *

del CalcYSE_HF

__version__ = get_versions()["version"]
del get_versions

__author__ = "Adrien Broquet"

__all__ = [
    "Conversion_Te_HF",
    "Curv_Moment",
    "Conversion_Tprofile_Te",
    "Planet_constants",
]
