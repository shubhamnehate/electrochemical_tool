# EIS models and equivalent circuit fitting

import numpy as np
from scipy.optimize import minimize, curve_fit
from scipy.interpolate import interp1d


# ── Circuit element impedances ─────────────────────────────────────────────────

def Z_R(R, omega):
    return np.full_like(omega, R, dtype=complex)

def Z_C(C, omega):
    return 1.0 / (1j * omega * C)

def Z_CPE(Q, n, omega):
    """Constant Phase Element: Z = 1 / (Q * (jω)^n).
    n=1 → ideal capacitor, n=0.5 → Warburg-like diffusion."""
    return 1.0 / (Q * (1j * omega) ** n)

def Z_Warburg(sigma, omega):
    """Semi-infinite Warburg diffusion: Z = σ(jω)^{-0.5}"""
    return sigma * (1j * omega) ** (-0.5)

def Z_Gerischer(Rg, tau, omega):
    """Gerischer element for porous electrode: Z = Rg / sqrt(1 + jω·τ)"""
    return Rg / np.sqrt(1 + 1j * omega * tau)


# ── Equivalent circuit models ──────────────────────────────────────────────────

def randles_impedance(omega, Rs, Rct, Cdl, sigma):
    """
    Randles circuit with semi-infinite Warburg diffusion.
    Topology: Rs + [Cdl || (Rct + Zw)]
    """
    Zw = Z_Warburg(sigma, omega)
    Zf = Rct + Zw
    Z_par = Zf / (1.0 + 1j * omega * Cdl * Zf)
    return Rs + Z_par


def randles_cpe_impedance(omega, Rs, Rct, Q, n, sigma):
    """
    Randles circuit with CPE replacing ideal capacitor.
    Topology: Rs + [CPE(Q,n) || (Rct + Zw)]
    More realistic for rough/porous electrode surfaces.
    """
    Zw = Z_Warburg(sigma, omega)
    Zf = Rct + Zw
    Z_cpe = Z_CPE(Q, n, omega)
    Y_par = 1.0 / Zf + 1.0 / Z_cpe
    return Rs + 1.0 / Y_par


def two_rc_impedance(omega, Rs, R1, C1, R2, C2):
    """
    Two-RC-loop circuit for SEI + charge transfer.
    Topology: Rs + [C1||R1] + [C2||R2]
    Used for battery electrodes with distinct SEI and double-layer processes.
    """
    Z1 = R1 / (1.0 + 1j * omega * R1 * C1)
    Z2 = R2 / (1.0 + 1j * omega * R2 * C2)
    return Rs + Z1 + Z2


def two_rc_warburg_impedance(omega, Rs, R1, C1, R2, C2, sigma):
    """
    Two-RC + Warburg for aged battery electrodes.
    SEI loop (R1,C1) + charge-transfer loop (R2,C2) + diffusion tail.
    """
    Z1 = R1 / (1.0 + 1j * omega * R1 * C1)
    Z2 = R2 / (1.0 + 1j * omega * R2 * C2)
    Zw = Z_Warburg(sigma, omega)
    return Rs + Z1 + Z2 + Zw


def pem_impedance(omega, Rmem, Rct, Cdl, sigma, L):
    """
    PEM fuel cell impedance model.
    Topology: jωL + Rmem + [Cdl || (Rct + Zw)]
    Includes membrane resistance, HOR/ORR kinetics, and mass-transport diffusion.
    """
    Zw = Z_Warburg(sigma, omega)
    Zf = Rct + Zw
    Z_par = Zf / (1.0 + 1j * omega * Cdl * Zf)
    Z_ind = 1j * omega * L
    return Z_ind + Rmem + Z_par


# ── Frequency generation ───────────────────────────────────────────────────────

def generate_frequencies(f_low=0.01, f_high=1e5, n_points=60):
    """Logarithmically spaced frequency array (Hz)."""
    return np.logspace(np.log10(f_low), np.log10(f_high), n_points)


# ── Degradation interpolation ──────────────────────────────────────────────────

def interpolate_params(fresh_params, aged_params, degradation_level):
    """
    Linearly interpolate circuit parameters between fresh (0) and fully aged (1).
    degradation_level: float in [0, 1]
    """
    d = np.clip(degradation_level, 0.0, 1.0)
    result = {}
    for key in fresh_params:
        result[key] = fresh_params[key] * (1 - d) + aged_params[key] * d
    return result


# ── Synthetic data generation with noise ──────────────────────────────────────

def generate_eis_data(model_fn, params, freqs, noise_level=0.02, seed=42):
    """
    Generate synthetic EIS spectra with Gaussian noise.
    noise_level: fraction of |Z| added as random noise.
    Returns complex impedance array.
    """
    rng = np.random.default_rng(seed)
    omega = 2 * np.pi * freqs
    Z_model = model_fn(omega, **params)
    noise_re = rng.normal(0, noise_level * np.abs(Z_model))
    noise_im = rng.normal(0, noise_level * np.abs(Z_model))
    return Z_model + (noise_re + 1j * noise_im)


# ── Equivalent circuit fitting ─────────────────────────────────────────────────

def fit_randles(freqs, Z_data, p0=None, bounds=None):
    """
    Fit Randles circuit to impedance data using nonlinear least squares.
    Returns fitted parameters dict and residual.
    """
    omega = 2 * np.pi * np.array(freqs)

    if p0 is None:
        p0 = [np.min(Z_data.real), 0.1, 1e-3, 0.02]
    if bounds is None:
        bounds = ([0, 1e-6, 1e-9, 1e-6], [10, 100, 10, 10])

    def model_flat(omega, Rs, Rct, Cdl, sigma):
        Z = randles_impedance(omega, Rs, Rct, Cdl, sigma)
        return np.concatenate([Z.real, Z.imag])

    Z_flat = np.concatenate([Z_data.real, Z_data.imag])
    try:
        popt, pcov = curve_fit(model_flat, omega, Z_flat, p0=p0, bounds=bounds,
                               maxfev=10000)
        Z_fit = randles_impedance(omega, *popt)
        residual = np.mean(np.abs(Z_data - Z_fit) / np.abs(Z_data)) * 100
        perr = np.sqrt(np.diag(pcov))
        return {
            "Rs": popt[0], "Rct": popt[1], "Cdl": popt[2], "sigma": popt[3],
            "Rs_err": perr[0], "Rct_err": perr[1],
            "Cdl_err": perr[2], "sigma_err": perr[3],
            "residual_pct": residual, "success": True,
            "Z_fit": Z_fit,
        }
    except Exception as e:
        return {"success": False, "error": str(e)}


def fit_two_rc_warburg(freqs, Z_data, p0=None):
    """Fit two-RC + Warburg circuit to impedance data."""
    omega = 2 * np.pi * np.array(freqs)

    if p0 is None:
        Rs_est = np.min(Z_data.real)
        p0 = [Rs_est, 0.02, 0.05, 0.08, 1e-3, 0.03]

    bounds = ([0]*6, [5, 5, 10, 5, 10, 10])

    def model_flat(omega, Rs, R1, C1, R2, C2, sigma):
        Z = two_rc_warburg_impedance(omega, Rs, R1, C1, R2, C2, sigma)
        return np.concatenate([Z.real, Z.imag])

    Z_flat = np.concatenate([Z_data.real, Z_data.imag])
    try:
        popt, pcov = curve_fit(model_flat, omega, Z_flat, p0=p0, bounds=bounds,
                               maxfev=20000)
        Z_fit = two_rc_warburg_impedance(omega, *popt)
        residual = np.mean(np.abs(Z_data - Z_fit) / np.abs(Z_data)) * 100
        perr = np.sqrt(np.diag(pcov))
        labels = ["Rs", "R_SEI", "C_SEI", "Rct", "Cdl", "sigma"]
        errs = [f"{l}_err" for l in labels]
        return {
            **dict(zip(labels, popt)),
            **dict(zip(errs, perr)),
            "residual_pct": residual, "success": True,
            "Z_fit": Z_fit,
        }
    except Exception as e:
        return {"success": False, "error": str(e)}


# ── DRT (Distribution of Relaxation Times) approximation ─────────────────────

def compute_drt_approx(freqs, Z_data, n_gaussian=5):
    """
    Simple DRT approximation: deconvolve the imaginary part of Z into
    a sum of Gaussians in log(τ) space.  Used to visualise overlapping
    time constants without a full Tikhonov regularisation.
    """
    omega = 2 * np.pi * np.array(freqs)
    tau = 1.0 / omega
    log_tau = np.log10(tau)
    Z_im = -Z_data.imag   # positive for capacitive arcs

    log_tau_range = np.linspace(log_tau.min(), log_tau.max(), 500)
    # Gaussian kernel density weighted by |Z_im|
    bw = 0.4
    drt = np.zeros(500)
    for lt, zi in zip(log_tau, Z_im):
        if zi > 0:
            drt += zi * np.exp(-0.5 * ((log_tau_range - lt) / bw) ** 2)
    drt /= (drt.max() + 1e-12)
    return log_tau_range, drt
