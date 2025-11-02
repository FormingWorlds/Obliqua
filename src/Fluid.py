import numpy as np
from scipy.fft import fft, fftshift
from scipy.interpolate import interp1d
from scipy.special import gamma


# -------------------- CONSTANTS --------------------
AU = 1.495978707e11  # m
G  = 6.6743e-11  # m^3 kg^-1 s^-2

# ---------------- FROM CONFIG PROTEUS ----------------
sigma_SB = 5.6704e-8  # W/m^2 K^4

# ---------------- FROM HELPFILE PROTEUS ----------------
    # Star
R_Sun = 6.95700e8  # m
M_Sun = 1.98847e30  # kg
L_Sun = 3.828e26  # W

S_mass   = 1 * M_Sun
S_radius = 1 * R_Sun
S_temp   = 5508.0  # K
S_surf_area = 4 * np.pi * (S_radius ** 2)
S_luminos   = S_surf_area * sigma_SB * (S_temp ** 4)
S_luminos2  = S_luminos / L_Sun

M_jup = 1898.13e24  # kg
R_jup = 69911e3  # m

M_Earth = 5.9724e24  # kg
R_Earth = 6371e3  # m

P_mass   = 8 * M_Earth
P_radius = 2 * R_Earth
P_grav_acc = G * P_mass / P_radius ** 2

rho_L = 2500.0  # kg/m^3 (lava)
rho_metal = 11e3  # kg/m^3 (core)

    # Orbit
spin_orbit_synchronized = 1

P_orb   = 1.0  # days
P_n_orb = 2.0 * np.pi / (P_orb * 3600.0 * 24.0)  # orbital freq (s^-1)
P_ecc   = 0.03

P_semimajoraxis = (G * S_mass / (P_n_orb ** 2)) ** (1.0 / 3.0)



# -------------------- MANTLE / PHYSICAL --------------------
mu = 3e11  # shear modulus Pa
Alpha = 0.3  # Andrade alpha
vis_solid = 1e21  # Pa s
tau_M = vis_solid / mu
tau_A = tau_M

emissivity = 0.9
albedo = 0.1

sigmaR = 10 ** (-3)
H_magma = P_radius / 50.0
R_solid_mantle = P_radius - H_magma
R_core = 0.4 * P_radius

volume_ocean = 4.0 / 3.0 * np.pi * (P_radius ** 3 - R_solid_mantle ** 3)
volume_core = 4.0 / 3.0 * np.pi * (R_core ** 3)
mass_ocean = volume_ocean * rho_L
mass_core = volume_core * rho_metal
mass_mantle = P_mass - mass_ocean - mass_core
rho_c = mass_mantle / (4.0 / 3.0 * np.pi * R_solid_mantle ** 3)
rho_ratio = rho_L / rho_c