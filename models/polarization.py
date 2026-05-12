# Polarization curve models — PEM fuel cells and Li-ion batteries

import numpy as np

R_GAS = 8.314    # J mol⁻¹ K⁻¹
F_CONST = 96485  # C mol⁻¹


# ── PEM Fuel Cell ──────────────────────────────────────────────────────────────

def nernst_voltage(T=353.15, P_H2=1.0, P_O2=0.21):
    """
    Reversible cell voltage from Nernst equation.
    T: temperature [K], P_H2/P_O2: partial pressures [atm]
    Returns E_rev [V].
    """
    E_std = 1.229 - 8.5e-4 * (T - 298.15)
    E_nernst = E_std + (R_GAS * T) / (2 * F_CONST) * np.log(P_H2 * np.sqrt(P_O2))
    return E_nernst


def activation_overpotential_pem(i, i0, T=353.15, alpha=0.5):
    """
    Activation overpotential using Tafel approximation (high-field Butler-Volmer).
    Valid for i >> i0.  Returns η_act [V].
    """
    i = np.maximum(i, 1e-9)
    return (R_GAS * T) / (alpha * F_CONST) * np.log(i / i0)


def ohmic_overpotential_pem(i, R_mem):
    """Ohmic loss: η_ohm = i · R_mem [V], R_mem in Ω·cm²"""
    return i * R_mem


def concentration_overpotential_pem(i, i_L, B=0.05):
    """
    Mass-transport (concentration) overpotential.
    B: empirical coefficient [V]; i_L: limiting current density [A cm⁻²].
    """
    ratio = np.clip(i / i_L, 0, 0.9999)
    return -B * np.log(1.0 - ratio)


def pem_polarization_curve(i_array, params):
    """
    Full PEM V-I polarization curve.
    params keys: T, P_H2, P_O2, i0, alpha, R_mem, i_L, B
    Returns dict with V_cell and individual overpotential arrays.
    """
    T = params.get("T", 353.15)
    P_H2 = params.get("P_H2", 1.0)
    P_O2 = params.get("P_O2", 0.21)
    i0 = params["i0"]
    alpha = params.get("alpha", 0.5)
    R_mem = params["R_mem"]
    i_L = params["i_L"]
    B = params.get("B", 0.05)

    E_rev = nernst_voltage(T, P_H2, P_O2)
    eta_act = activation_overpotential_pem(i_array, i0, T, alpha)
    eta_ohm = ohmic_overpotential_pem(i_array, R_mem)
    eta_conc = concentration_overpotential_pem(i_array, i_L, B)
    V_cell = E_rev - eta_act - eta_ohm - eta_conc

    # Clamp to physical range
    V_cell = np.clip(V_cell, 0, None)

    return {
        "V_cell": V_cell,
        "E_rev": np.full_like(i_array, E_rev),
        "eta_act": eta_act,
        "eta_ohm": eta_ohm,
        "eta_conc": eta_conc,
        "P_density": i_array * V_cell,   # W cm⁻²
    }


def pem_degradation_params(base_params, deg_level):
    """
    Return degraded PEM parameters at degradation_level ∈ [0,1].
    Models: catalyst dissolution (↓i0), membrane thinning (↑R_mem),
    flooding/drying (↓i_L).
    """
    d = np.clip(deg_level, 0, 1)
    p = dict(base_params)
    p["i0"]   = base_params["i0"]   * (1 - 0.85 * d)   # ECSA loss
    p["R_mem"] = base_params["R_mem"] * (1 + 2.0 * d)   # membrane resistance increase
    p["i_L"]  = base_params["i_L"]  * (1 - 0.40 * d)   # mass-transport degradation
    p["alpha"] = base_params["alpha"] * (1 - 0.10 * d)  # kinetic parameter shift
    return p


# ── Li-ion Battery ─────────────────────────────────────────────────────────────

# OCV models for common electrode materials

def ocv_nmc(soc):
    """NMC811 cathode OCV vs. Li/Li+ [V] — polynomial fit to Ruan et al. 2021."""
    s = np.clip(soc, 0.01, 0.99)
    return (3.85
            + 0.60 * s
            - 1.20 * s**2
            + 0.80 * s**3
            - 0.28 * np.log(s)
            + 0.18 * np.log(1 - s))


def ocv_lfp(soc):
    """LFP cathode OCV [V] — plateau model (Dreyer et al.)."""
    s = np.clip(soc, 0.02, 0.98)
    # Sigmoid transitions around the two-phase plateau at ~3.45 V
    V_plateau = 3.42
    delta = 0.015
    V = V_plateau + delta * np.tanh((s - 0.5) / 0.12)
    # Slope at ends
    V += 0.15 * (s - 0.5)**3
    return V


def ocv_graphite(soc):
    """Graphite anode OCV vs. Li/Li+ [V] — multi-step staging model."""
    s = np.clip(soc, 0.01, 0.99)
    # Approximate staging plateaus
    V = (0.07 + 0.07 * np.exp(-70 * s)
         + 0.17 * np.exp(-5.0 * s)
         + 0.10 * s
         - 0.045 * np.log(s)
         + 0.02 * np.log(1 - s))
    return V


def battery_voltage(soc_array, params, direction="discharge"):
    """
    Li-ion cell voltage profile during (dis)charge.
    params keys: Q_cell, C_rate, R_int, cathode, T
    Returns dict with V, I, P, individual contributions.
    """
    cathode = params.get("cathode", "NMC811")
    C_rate = params.get("C_rate", 1.0)
    Q_cell = params.get("Q_cell", 200.0)   # mAh/g
    R_int = params.get("R_int", 0.020)     # Ω (normalised to cell)
    T = params.get("T", 298.15)

    # Current [A/g], sign convention: discharge positive
    I = C_rate * Q_cell / 1000.0   # A/g (scaled)

    # OCV
    ocv_fn = {"NMC811": ocv_nmc, "LFP": ocv_lfp, "Graphite": ocv_graphite}.get(cathode, ocv_nmc)
    V_ocv = ocv_fn(soc_array)

    # Activation overpotential (Butler-Volmer linearised for small η)
    i0_cell = params.get("i0", 5.0)    # mA/g
    alpha = params.get("alpha", 0.5)
    F_norm = F_CONST / 1000.0          # keep units consistent
    eta_act = (R_GAS * T / (alpha * F_CONST)) * np.arcsinh(I / (2 * i0_cell * 1e-3))

    # Ohmic drop
    eta_ohm = I * R_int

    if direction == "discharge":
        V_cell = V_ocv - eta_act - eta_ohm
    else:
        V_cell = V_ocv + eta_act + eta_ohm

    V_cell = np.clip(V_cell, 2.5, 4.5)

    return {
        "V_cell": V_cell,
        "V_ocv": V_ocv,
        "eta_act": np.full_like(soc_array, eta_act),
        "eta_ohm": np.full_like(soc_array, eta_ohm),
        "P_cell": V_cell * I,
        "I": I,
    }


def battery_degradation_params(base_params, deg_level):
    """
    Return degraded battery parameters at degradation_level ∈ [0,1].
    Models: Li plating/SEI (↑R_int), active material loss (↓i0),
    lithium inventory loss (SOC window shrinks).
    """
    d = np.clip(deg_level, 0, 1)
    p = dict(base_params)
    p["R_int"] = base_params["R_int"] * (1 + 4.0 * d)    # SEI resistance growth
    p["i0"]    = base_params["i0"]    * (1 - 0.60 * d)   # active site loss
    p["Q_cell"] = base_params["Q_cell"] * (1 - 0.35 * d) # capacity fade
    return p


# ── dQ/dV analysis ─────────────────────────────────────────────────────────────

def compute_dqdv(soc_array, V_array, Q_total):
    """
    Compute differential capacity dQ/dV from voltage-SOC curve.
    Returns voltage array and dQ/dV array [Ah V⁻¹].
    """
    Q_array = soc_array * Q_total
    dQ = np.diff(Q_array)
    dV = np.diff(V_array)
    # Avoid division by near-zero dV
    mask = np.abs(dV) > 1e-6
    dqdv = np.zeros(len(dV))
    dqdv[mask] = np.abs(dQ[mask] / dV[mask])
    V_mid = 0.5 * (V_array[:-1] + V_array[1:])
    return V_mid, dqdv
