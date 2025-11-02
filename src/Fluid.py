import numpy as np
from scipy.fft import fft, fftshift
from scipy.interpolate import interp1d



# -------------------- CONSTANTS --------------------
AU = 1.495978707e11  # m
G  = 6.6743e-11  # m^3 kg^-1 s^-2

# ---------------- Degree Love number ---------------
n = 2
m = 2
k_min = -30
k_max = 40
k_range = np.arange(k_min, k_max + 1)

# ---------------- FROM HELPFILE PROTEUS ----------------
    # Star
M_Sun = 1.98847e30  # kg
S_mass   = 1 * M_Sun

M_Earth    = 5.9724e24  # kg
P_mass     = 8 * M_Earth

sigmaR = 10 ** (-3)         # friction between solid and liquid, should prob be varied

def get_fluid_k2(omega, axial, ecc, sma, density, radius, visc, shear, bulk, visc_thresh, N_sigma=301):
    
    P_radius   = np.max(radius)
    P_grav_acc = G * P_mass / P_radius ** 2

    r_l, r_s, rho_l, rho_s = find_liquid_region_and_densities(density, radius, visc, visc_thresh)
    H_magma = P_radius - np.max(r_s)
    
    rho_ratio = rho_l / rho_s

    P_n_orb = omega # orbital freq (s^-1)
    Omega   = axial # axial freq (s^-1)

    # get hansen coefficients
    k_range2, X_hansen = get_hansen(ecc)

    # set dissipative frequency range
    T_range = np.logspace(-15, 6, N_sigma)  # periods
    sigma_range = 2.0 * np.pi / (T_range * 1e3 * 365.25 * 24.0 * 3600.0)
    sigma_range = sigma_range.reshape(-1)

    # preallocate (complex for viscoelastic)
    k_T_homo = np.zeros((n, N_sigma), dtype=complex)
    k_L_homo = np.zeros((n, N_sigma), dtype=complex)

    # could include Andrade solid tides here, instead of propagator method
    # ...

    # get 2,2 harmonic
    k_T_22_homo = k_T_homo[1, :].reshape(-1)
    k_L_22_homo = k_L_homo[1, :].reshape(-1)

    # fluid Love Numbers
    k22_fluid_high_friction, k22_total = compute_fluid_lovenumbers(
        n,
        sigma_range,
        k_T_22_homo,
        k_L_22_homo,
        rho_ratio,
        P_grav_acc,         # !!!
        H_magma,            # !!!
        sigmaR,
        P_radius,           # !!!
    )

    # interpolate k2 love number arrays
    mu_n = n * (n + 1)
    ksi_n = 3.0 / (2.0 * n + 1.0) * rho_ratio
    sigP_n = np.sqrt(mu_n * P_grav_acc * H_magma / P_radius ** 2)

    for kk in range(N_sigma):
        sigma = sigma_range[kk]
        sigT = sigma - 1j * sigmaR
        k22_fluid_high_friction[kk] = -ksi_n * sigP_n ** 2 / (sigma * sigT - sigP_n ** 2)
        k22_total[kk] = k_T_22_homo[kk] + (1.0 + k_L_22_homo[kk]) * k22_fluid_high_friction[kk]

    k22_total = k22_total.reshape(-1)

    # Build symmetric full spectrum for interpolation (like MATLAB concatenation)
    full_sigma_range = np.concatenate((-sigma_range, np.flip(sigma_range)))
    full_k22_total = np.concatenate((-k22_total, np.flip(k22_total)))
    full_k22_homo = np.concatenate((-k_T_22_homo, np.flip(k_T_22_homo)))

    imag_full_k22 = np.imag(full_k22_total)
    imag_solid_k22 = np.imag(full_k22_homo)

    # interpolation functions for imaginary parts (extrapolate outside)
    imag_k22_full_spectrum = interp1d(full_sigma_range, imag_full_k22, kind="linear", fill_value="extrapolate")
    imag_k22_solid_spectrum = interp1d(full_sigma_range, imag_solid_k22, kind="linear", fill_value="extrapolate")

    # calculate tidal heating
    A_22k_e = np.zeros(len(k_range), dtype=float)
    U_22k_e = np.zeros(len(k_range), dtype=complex)
    P_T_k_total = np.zeros(len(k_range), dtype=float)
    P_T_k_solid = np.zeros(len(k_range), dtype=float)

    for ikk, kk in enumerate(k_range):
        sigma = 2.0 * Omega - kk * P_n_orb
        A_22k_e[ikk] = np.sqrt(6.0 * np.pi / 5.0) * X_hansen[ikk]  # Eq.33
        U_22k_e[ikk] = (G * S_mass / sma) * (P_radius / sma) ** 2 * A_22k_e[ikk]  # Eq.32

        img_full_k22 = float(imag_k22_full_spectrum(sigma))
        img_solid_k22 = float(imag_k22_solid_spectrum(sigma))

        prefactor = 5.0 * P_radius * sigma / (8.0 * np.pi * G)
        U2 = abs(U_22k_e[ikk]) ** 2

        P_T_k_total[ikk] = prefactor * img_full_k22 * U2
        P_T_k_solid[ikk] = prefactor * img_solid_k22 * U2

    P_tidal_total = -np.sum(P_T_k_total)
    P_tidal_solid = -np.sum(P_T_k_solid)

    return P_T_k_total, P_tidal_total, img_full_k22


# -------------------- Helper functions --------------------

def find_liquid_region_and_densities(density, radius, visc, visc_thresh):
    """
    Identify liquid shells (viscosity < visc_thresh) and compute
    mass-weighted mean density of the liquid and solid regions.

    Parameters
    ----------
    density : array_like
        Density of each radial shell (kg/m^3).
    radius : array_like
        Outer radius of each shell (m). Must be increasing.
    visc : array_like
        Dynamic viscosity of each shell (Pa s).
    visc_thresh : float
        Threshold viscosity for liquid behaviour.

    Returns
    -------
    r_liquid : ndarray
        Radii of shells where visc < visc_thresh.
    r_solid : ndarray
        Radii of shells where visc >= visc_thresh.
    mean_rho_liquid : float
        Mass-weighted mean density in liquid region. NaN if region empty.
    mean_rho_solid : float
        Mass-weighted mean density in solid region. NaN if region empty.
    """

    density = np.asarray(density, dtype=float)
    radius  = np.asarray(radius, dtype=float)
    visc    = np.asarray(visc, dtype=float)

    # compute shell volumes & masses
        # interpret radius[i] as the OUTER boundary of shell i.
        # the inner boundary is radius[i-1], with r_inner = 0 for i=0.
    r_outer = radius
    r_inner = np.zeros_like(radius)
    r_inner[1:] = radius[:-1]    # previous radius

    shell_vol = (4.0/3.0) * np.pi * (r_outer**3 - r_inner**3)
    shell_mass = density * shell_vol

    # Define liquid region mask
    mask_liquid = visc < visc_thresh
    mask_solid = ~mask_liquid

    r_liquid = radius[mask_liquid]
    r_solid  = radius[mask_solid]

    # mass-weighted mean densities
    def mass_weighted_mean_rho(mask):
        if not np.any(mask):
            return np.nan
        m = shell_mass[mask]
        return np.sum(density[mask] * m) / np.sum(m)

    mean_rho_liquid = mass_weighted_mean_rho(mask_liquid)
    mean_rho_solid  = mass_weighted_mean_rho(mask_solid)

    return r_liquid, r_solid, mean_rho_liquid, mean_rho_solid


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

def get_hansen(ecc):   

    k_range2, X_hansen = hansen_fft(-(n + 1), m, ecc, k_min, k_max, 2 ** 18)
    X_hansen = X_hansen.reshape(-1)  # 1D array

    return k_range2, X_hansen


# -------------------- FLUID TIDE K2 LOVE NUMBERS --------------------

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

