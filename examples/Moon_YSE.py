from Te_HF_Conversion import Conversion_Te_HF, Planet_constants
import numpy as np
from astropy.constants import G as _G

GM = 4902.80007e9
R_mpr = 1737.151e3
mass = GM / _G.value
g = GM / R_mpr**2
rhobar = 3 * mass / (np.pi * 4 * R_mpr**3)
# From Mohit & Phillips
Q_olivine = 540e3
A_olivine = 2.4e-16
n_olivine = 3.5
Q_diabase = 485e3
A_diabase = 1.2e-26
n_diabase = 4.7
R_gas = Planet_constants.R_gas

## Arbitrary constants
rhoc = 2550.0  # density of the crust (kg m-3)
rhom = 3220.0  # density of the mantle (kg m-3)
Ts = 250.0  # surface temperature (K)
young = 1e11  # young's modulus (Pa)
poisson = 0.25  # poisson's ratio
k_crust = 1.5  # thermal conductivity of the crust (W m−1 K−1)
k_mantle = 3.0  # thermal conductivity of the mantle (W m−1 K−1)
sig_y = 10e6  # Bounding stress (Pa)
H_c = 0.0  # 4.9e-11  # Average volumetric crustal heat production (W kg-1)
Te = 90e3  # elastic thickness (m)
Tc = 40e3  # crustal thickness (m)
K_curv = 4.1e-10  # plate curvature (m-1)
eps = 1e-17  # strain rate (s-1)

# Heat flow bounds to be investigated to improve speed convergence
# (if None, automatically bounded, based on the input and expected results)
HF_min = None
HF_max = None
plot_YSE = False  # Plot yield strength envelope on the fly for each tested heat flow
plot = True  # Plot the final yield strength envelope and associated temperature profile

# More input parameters can be input to increase speed (see doc)
kwargs = dict(quiet=False, plot=plot, plot_YSE=plot_YSE, HF_min=HF_min, HF_max=HF_max)

Conversion_Te_HF(
    Q_diabase,
    A_diabase,
    n_diabase,
    Q_olivine,
    A_olivine,
    n_olivine,
    g,
    rhobar,
    R_mpr,
    Te,
    Tc,
    k_crust,
    k_mantle,
    K_curv,
    H_c,
    rhoc,
    rhom,
    eps,
    young,
    poisson,
    sig_y,
    Ts,
    R_gas,
    **kwargs,
)
