import numpy as np

# gas cst CODATA (J mol-1 K-1)
R_gas = 8.314462618
# grav cst CODATA (m3 kg-1 s-2)
G_cst = 6.6743e-11

# Rheological constants from Grott & Breuer 2008

# activation nrj dry diabase (J mol-1)
Q_diabase_dry = 4.88e5
# empirical cst dry diabase (Pa-n s-1)
A_diabase_dry = 1.1e-26
# rheological cst dry diabase
n_diabase_dry = 4.7

# activation nrj wet diabase (J mol-1)
Q_diabase_wet = 2.76e5
# empirical cst wet diabase (Pa-n s-1)
A_diabase_wet = 3.1e-20
# rheological cst wet diabase
n_diabase_wet = 3.05

# activation nrj dry olivine (J mol-1)
Q_olivine_dry = 5.4e5
# empirical cst dry olivine (Pa-n s-1)
A_olivine_dry = 2.4e-16
# rheological cst dry olivine
n_olivine_dry = 3.5

# activation nrj wet olivine
Q_olivine_wet = 4.2e5
# empirical cst wet olivine
A_olivine_wet = 1.9e-15
# rheological cst wet olivine
n_olivine_wet = 3.0

# Mars constants
# gm (m3 s-2) – Konopliv et al. (2016)
Mars_gm = 0.4282837581575610e14
# Mean planetary radius (m) – MarsTopo2600: Wieczorek, M. A. (2015)
Mars_mpr = 3389.500e3
# Mass (kg)
Mars_mass = Mars_gm / G_cst
# Gravitational attraction (m s-2)
Mars_g0 = Mars_gm / Mars_mpr**2
# Average density  (kg m-3)
Mars_density = Mars_mass * 3.0 / 4.0 / np.pi / (Mars_mpr**3)

# Venus constants
# gm (m3 s-2) – Konopliv et al. (1999)
Venus_gm = 324858592079000.0
# Mean planetary radius (m) – VenusTopo719: Wieczorek, M. A. (2015)
Venus_mpr = 6051.878e3
# Mass (kg)
Venus_mass = Venus_gm / G_cst
# Gravitational attraction (m s-2)
Venus_g0 = Venus_gm / Venus_mpr**2
# Average density (kg m-3)
Venus_density = Venus_mass * 3.0 / 4.0 / np.pi / (Venus_mpr**3)
