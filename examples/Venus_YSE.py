from Te_HF_Conversion import *

plot_YSE = True  # Plot yield strength envelope on the fly for each tested heat flow
plot = True  # Plot the final yield strength envelope and associated temperature profile

g = Planet_constants.Venus_g0
rhobar = Planet_constants.Venus_density
R_mpr = Planet_constants.Venus_mpr
Q_olivine = Planet_constants.Q_olivine_dry
A_olivine = Planet_constants.A_olivine_dry
n_olivine = Planet_constants.n_olivine_dry
Q_diabase = Planet_constants.Q_diabase_dry
A_diabase = Planet_constants.A_diabase_dry
n_diabase = Planet_constants.n_diabase_dry
R_gas = Planet_constants.R_gas

## Arbitrary constants
rhoc = 2800.0  # density of the crust (kg m-3)
rhom = 3300.0  # density of the mantle (kg m-3)
Ts = 740.0  # surface temperature (K)
young = 1e11  # young's modulus (Pa)
poisson = 0.25  # poisson's ratio
k_crust = 3.0  # thermal conductivity of the crust (W m−1 K−1)
k_mantle = 4.0  # thermal conductivity of the mantle (W m−1 K−1)
sig_y = 50e6  # Bounding stress (Pa)
H_c = 0.  # Average volumetric crustal heat production (W kg-1)
Te = 100e3  # elastic thickness (m)
Tc = 70e3  # crustal thickness (m)
K_curv = 1e-7  # plate curvature (m-1)
eps = 1e-16  # strain rate (s-1)

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
    quiet=False,
    plot=plot,
    plot_YSE=plot_YSE,
    HF_min=5,
    HF_max=15,
)
