import numpy as np
import matplotlib.pyplot as plt
from Te_HF_Conversion import *

g = Planet_constants.Mars_g0
rhobar = Planet_constants.Mars_density
R_mpr = Planet_constants.Mars_mpr
Q_olivine = Planet_constants.Q_olivine_wet
A_olivine = Planet_constants.A_olivine_wet
n_olivine = Planet_constants.n_olivine_wet
Q_diabase = Planet_constants.Q_diabase_wet
A_diabase = Planet_constants.A_diabase_wet
n_diabase = Planet_constants.n_diabase_wet
R_gas = Planet_constants.R_gas

## Arbitrary constants
rhoc = 2800.0  # density of the crust (kg m-3)
rhom = 3500.0  # density of the mantle (kg m-3)
young = 1e11  # young's modulus (Pa)
poisson = 0.25  # poisson's ratio
sig_y = 50e6  # Bounding stress (Pa)
Tc = 50e3  # crustal thickness (m)
eps = 1e-16  # strain rate (s-1)

max_depth = 200  # Maximum depth investigated (km)
step_depth = 0.1  # Depth steps (km)
thermal_gradient = 10  # Thermal gradient (K / km)
z_profile = np.arange(0, max_depth, step=step_depth) * 1e3  # to meters
T_profile = 210 + z_profile * thermal_gradient / 1e3

# An example temperature profile with a thermal anomaly
depth_anomaly = 20  # depth of the thermal anomaly (km)
thickness_anomaly = 30  # Thickness of the thermal anomaly (km)
max_incre = int(
    np.min([thickness_anomaly - 1, max_depth - thickness_anomaly]) / step_depth
)
dT = 350  # Temperature contrast of the thermal anomaly
T_profile[
    int(depth_anomaly / step_depth) : int(depth_anomaly / step_depth) + max_incre
] += (np.cos(np.linspace(-np.pi / 2.0, np.pi / 2.0, max_incre)) * dT)

# Elastic thickness (km) bounds to be investigated to improve speed convergence
# (defaults are 5 to 300 km)
Te_min = 20
Te_max = 40
plot = True  # Plot the final yield strength envelope and associated temperature profile

# More input parameters can be input to increase speed (see doc)
kwargs = dict(quiet=False, plot=plot, Te_min=Te_min, Te_max=Te_max)

Conversion_Tprofile_Te(
    T_profile,
    z_profile,
    Q_diabase,
    A_diabase,
    n_diabase,
    Q_olivine,
    A_olivine,
    n_olivine,
    g,
    rhobar,
    R_mpr,
    Tc,
    rhoc,
    rhom,
    eps,
    young,
    poisson,
    sig_y,
    R_gas,
    **kwargs,
)
