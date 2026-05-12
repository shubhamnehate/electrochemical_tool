# Capacity fade and ageing dynamics models

import numpy as np

R_GAS = 8.314   # J mol⁻¹ K⁻¹


# ── Calendar ageing models ─────────────────────────────────────────────────────

def calendar_fade_sqrt(t_days, Q0, k_cal, beta=0.5):
    """
    Square-root (diffusion-controlled SEI growth) calendar capacity fade.
    Q(t) = Q0 * (1 - k_cal * t^beta)
    Validated for NMC/graphite at elevated temperatures (Keil et al. 2016).
    """
    return Q0 * (1.0 - k_cal * t_days ** beta)


def calendar_fade_arrhenius(t_days, Q0, A_cal, Ea, T_K, SOC=0.5):
    """
    Arrhenius-weighted calendar fade.  k = A * exp(-Ea/RT) * f(SOC)
    f(SOC) = exp(0.5*SOC) models higher fade at elevated state of charge.
    """
    k_T = A_cal * np.exp(-Ea / (R_GAS * T_K))
    soc_factor = np.exp(0.5 * SOC)
    return Q0 * (1.0 - k_T * soc_factor * t_days ** 0.5)


# ── Cycle ageing models ────────────────────────────────────────────────────────

def cycle_fade_linear(n_cycles, Q0, k_cyc):
    """Linear cycle fade: Q(n) = Q0 * (1 - k_cyc * n)"""
    return np.maximum(Q0 * (1.0 - k_cyc * n_cycles), 0.0)


def cycle_fade_power(n_cycles, Q0, k_cyc, beta=0.55):
    """
    Power-law cycle fade: Q(n) = Q0 * (1 - k_cyc * n^beta)
    beta < 1 captures decelerating fade typical of NMC cells.
    """
    return np.maximum(Q0 * (1.0 - k_cyc * n_cycles ** beta), 0.0)


def cycle_fade_exponential(n_cycles, Q0, k_cyc):
    """Exponential cycle fade: Q(n) = Q0 * exp(-k_cyc * n)"""
    return Q0 * np.exp(-k_cyc * n_cycles)


def cycle_fade_two_stage(n_cycles, Q0, k1, k2, n_knee):
    """
    Two-stage (knee-point) fade model.
    Stage 1 (n < n_knee): slow linear fade
    Stage 2 (n >= n_knee): accelerated exponential fade (Li plating onset).
    """
    Q = np.zeros_like(n_cycles, dtype=float)
    mask1 = n_cycles < n_knee
    mask2 = ~mask1
    Q[mask1] = Q0 * (1.0 - k1 * n_cycles[mask1])
    Q_knee = Q0 * (1.0 - k1 * n_knee)
    Q[mask2] = Q_knee * np.exp(-k2 * (n_cycles[mask2] - n_knee))
    return np.maximum(Q, 0.0)


# ── Combined ageing (Wang et al. model) ───────────────────────────────────────

def combined_fade_wang(n_cycles, t_days, Q0, B, Ea, T_K, C_rate, z=0.55):
    """
    Wang et al. (2011) combined capacity fade model.
    Q_loss = B * exp(-Ea/RT) * (Ah_throughput)^z
    Ah_throughput = n_cycles * C_rate * Q0
    Used widely for LFP and NMC cells.
    """
    Ah = n_cycles * C_rate * Q0 / 1000.0   # convert mAh/g → Ah
    k = B * np.exp(-Ea / (R_GAS * T_K))
    Q_loss = k * Ah ** z
    return np.maximum(Q0 - Q_loss, 0.0)


# ── PEM fuel cell performance degradation ─────────────────────────────────────

def pem_voltage_decay(operating_hours, V0, k_deg, k_startup=0.0, n_startups=0):
    """
    PEM fuel cell voltage decay over lifetime.
    k_deg: steady-state degradation [μV/h] — typically 3–10 μV/h (DOE target: <1 μV/h)
    k_startup: voltage loss per startup [μV/startup] — membrane pinhole growth.
    """
    V_decay = k_deg * 1e-6 * operating_hours + k_startup * 1e-6 * n_startups
    return np.maximum(V0 - V_decay, 0.0)


def pem_ecsa_loss(n_cycles, ECSA0, k_diss, k_osten, P_O2=0.21):
    """
    Pt catalyst ECSA (electrochemically active surface area) loss.
    Combined Pt dissolution + Ostwald ripening model.
    k_diss: dissolution rate constant [cm² / cycle]
    k_osten: Ostwald ripening coefficient
    """
    ECSA = ECSA0 * np.exp(-(k_diss + k_osten) * n_cycles * P_O2 ** 0.5)
    return np.maximum(ECSA, 0.0)


# ── Temperature & C-rate sensitivity ──────────────────────────────────────────

def temperature_acceleration_factor(T_K, T_ref_K=298.15, Ea=25000):
    """
    Arrhenius acceleration factor relative to reference temperature.
    Higher T → faster ageing.  Ea ≈ 25 kJ/mol (typical for NMC SEI growth).
    """
    return np.exp((Ea / R_GAS) * (1.0 / T_ref_K - 1.0 / T_K))


def crate_acceleration_factor(C_rate, C_ref=1.0, gamma=1.3):
    """
    C-rate acceleration factor for cycle ageing.
    Higher C-rate → more mechanical stress and heat → faster fade.
    gamma: empirical exponent (typically 1.2–1.5 for NMC).
    """
    return (C_rate / C_ref) ** gamma


def effective_fade_rate(base_k, T_K, C_rate,
                        T_ref_K=298.15, Ea=25000, gamma=1.3):
    """Combine temperature and C-rate acceleration on a base fade rate."""
    af_T = temperature_acceleration_factor(T_K, T_ref_K, Ea)
    af_C = crate_acceleration_factor(C_rate)
    return base_k * af_T * af_C


# ── Resistance growth models ───────────────────────────────────────────────────

def sei_resistance_growth(n_cycles, R0, k_sei, beta=0.5):
    """
    SEI layer resistance growth: R(n) = R0 * (1 + k_sei * n^beta)
    Follows square-root kinetics (diffusion-limited SEI growth).
    """
    return R0 * (1.0 + k_sei * n_cycles ** beta)


def ohmic_resistance_growth_pem(hours, R0, k_mem, k_corr):
    """
    PEM membrane + contact resistance growth over time.
    k_mem: membrane hydration-related drift [Ω·cm²/h]
    k_corr: corrosion-driven contact resistance increase [Ω·cm²/h^0.5]
    """
    return R0 + k_mem * hours + k_corr * np.sqrt(hours)
