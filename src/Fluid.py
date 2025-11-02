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


# -------------------- Helper functions --------------------

def nextpow2(x):
    """Return the exponent p such that 2**p >= x (integer p)."""
    return int(np.ceil(np.log2(x))) if x > 0 else 0


def kepler_newton(M, e):
    """
    Solve Kepler's equation E - e*sin(E) = M for eccentric anomaly E
    using Newton iterations. 
    """
    M = np.asarray(M, dtype=float)
    E = M.copy()
    if e > 0:
        E = M + e * np.sin(M) / (1 - np.sin(M + e) + np.sin(M))  # Danby-like tweak

    for _ in range(10):
        f = E - e * np.sin(E) - M
        fp = 1 - e * np.cos(E)
        dE = -f / fp
        E = E + dE
        if np.max(np.abs(dE)) < 1e-13:
            break

    return np.mod(E, 2 * np.pi)

def hansen_fft(n, m, e, kmin, kmax, N=None):
    """
    Hansen coefficients X_k^{n,m}(e) via FFT on mean anomaly.
    Returns (k, Xkm) where Xkm is real-valued (imag ~ roundoff).
    """
    # Adaptive N selection (power of 2)
    if N is None:
        width = max(64, 4 * (kmax - kmin + 1))
        target = width * max(8, int(np.ceil(16 / (1 - e + np.finfo(float).eps))))
        p = max(12, int(np.ceil(np.log2(target))))
        N = 2 ** p
    else:
        p = nextpow2(N)
        if 2 ** p != N:
            N = 2 ** p

    # Mean anomaly grid
    M = np.arange(N) * (2 * np.pi / N)

    # Solve Kepler
    E = kepler_newton(M, e)

    ce = np.cos(E)
    se = np.sin(E)
    r_over_a = 1 - e * ce
    v = np.arctan2(np.sqrt(1 - e ** 2) * se, ce - e)

    f = (r_over_a ** n) * np.exp(1j * m * v)

    F = fft(f) / N
    F = fftshift(F)

    k_all = np.arange(-N // 2, N // 2)
    mask = (k_all >= kmin) & (k_all <= kmax)
    k = k_all[mask]
    Zk = F[mask]
    Xkm = np.real(Zk)

    return k, Xkm

# -------------------- HANSEN COEFFICIENTS --------------------
n = 2
m = 2
k_min = -30
k_max = 40
k_range = np.arange(k_min, k_max + 1)

k_range2, X_hansen = hansen_fft(-(n + 1), m, P_ecc, k_min, k_max, 2 ** 18)
X_hansen = X_hansen.reshape(-1)  # 1D array


# -------------------- FLUID TIDE K2 LOVE NUMBER --------------------

def compute_fluid_lovenumbers(
    n,
    sigma_range,
    k_T_22_homo,
    k_L_22_homo,
    rho_ratio,
    P_grav_acc,
    H_magma,
    sigmaR,
    P_radius,
):
    """
    Compute magma-ocean fluid Love numbers for degree n=2 using Eq(7) and Eq(28).

    Returns:
        k22_fluid_high_friction (N_sigma,)
        k22_total (N_sigma,)
    """
    N_sigma = len(sigma_range)

    mu_n = n * (n + 1)
    ksi_n = 3 / (2 * n + 1) * rho_ratio
    sigP_n = np.sqrt(mu_n * P_grav_acc * H_magma / P_radius**2)

    k22_fluid = np.zeros(N_sigma, dtype=complex)
    k22_total = np.zeros(N_sigma, dtype=complex)

    for i, sigma in enumerate(sigma_range):
        sigT = sigma - 1j * sigmaR

        # Eq. 7
        k22_fluid[i] = -ksi_n * sigP_n**2 / (sigma * sigT - sigP_n**2)

        # Eq. 28 (ignoring crust terms)
        k22_total[i] = k_T_22_homo[i] + (1 + k_L_22_homo[i]) * k22_fluid[i]

    return k22_fluid, k22_total

    # set dissipative frequency range
N_sigma = 301
T_range = np.logspace(-15, 6, N_sigma)  # periods
sigma_range = 2.0 * np.pi / (T_range * 1e3 * 365.25 * 24.0 * 3600.0)
sigma_range = sigma_range.reshape(-1)

    # preallocate (complex for viscoelastic)
k_T_homo = np.zeros((n, N_sigma), dtype=complex)
k_L_homo = np.zeros((n, N_sigma), dtype=complex)

    # get 2,2 harmonic
k_T_22_homo = k_T_homo[1, :].reshape(-1)
k_L_22_homo = k_L_homo[1, :].reshape(-1)


    # --- Fluid Love Numbers ---
k22_fluid_high_friction, k22_total = compute_fluid_lovenumbers(
    n,
    sigma_range,
    k_T_22_homo,
    k_L_22_homo,
    rho_ratio,
    P_grav_acc,
    H_magma,
    sigmaR,
    P_radius,
)
