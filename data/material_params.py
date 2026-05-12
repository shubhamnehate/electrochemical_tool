# Electrochemical parameter database — PEM and Li-ion electrode materials

# ─────────────────────────────────────────────────────────────────────────────
# EIS equivalent circuit parameters (Randles circuit)
# Rs [Ω]  Rct [Ω]  Cdl [F]  sigma [Ω·s^-0.5]
# ─────────────────────────────────────────────────────────────────────────────

EIS_PARAMS = {
    "NMC811 Cathode": {
        "circuit": "randles",
        "fresh":   {"Rs": 0.022, "Rct": 0.075, "Cdl": 0.48, "sigma": 0.022},
        "aged":    {"Rs": 0.065, "Rct": 0.280, "Cdl": 0.18, "sigma": 0.090},
        "unit_area": "Ω (normalised to electrode area)",
        "notes": "NMC811 // graphite full cell, 25°C, C/10 SOC=50%",
    },
    "LFP Cathode": {
        "circuit": "randles",
        "fresh":   {"Rs": 0.012, "Rct": 0.035, "Cdl": 0.82, "sigma": 0.015},
        "aged":    {"Rs": 0.028, "Rct": 0.095, "Cdl": 0.60, "sigma": 0.040},
        "unit_area": "Ω (normalised to electrode area)",
        "notes": "LFP // graphite 26650 cell, 25°C, SOC=50%",
    },
    "Graphite Anode": {
        "circuit": "two_rc_warburg",
        "fresh":   {
            "Rs": 0.015, "R_SEI": 0.018, "C_SEI": 0.12,
            "Rct": 0.055, "Cdl": 0.32, "sigma": 0.018
        },
        "aged":    {
            "Rs": 0.025, "R_SEI": 0.085, "C_SEI": 0.08,
            "Rct": 0.160, "Cdl": 0.18, "sigma": 0.055
        },
        "unit_area": "Ω (normalised to electrode area)",
        "notes": "Graphite half-cell vs Li/Li+, 25°C, 1M LiPF6 EC/DMC",
    },
    "Pt/C PEM Cathode": {
        "circuit": "pem",
        "fresh":   {"Rs": 0.00, "Rmem": 0.14, "Rct": 0.11, "Cdl": 0.035, "sigma": 0.08, "L": 3e-7},
        "aged":    {"Rs": 0.00, "Rmem": 0.28, "Rct": 0.38, "Cdl": 0.015, "sigma": 0.20, "L": 3e-7},
        "unit_area": "Ω·cm²",
        "notes": "PEM single cell, H2/air, 80°C, 100% RH, 1 bar",
    },
    "NCA Cathode": {
        "circuit": "randles",
        "fresh":   {"Rs": 0.018, "Rct": 0.062, "Cdl": 0.55, "sigma": 0.020},
        "aged":    {"Rs": 0.058, "Rct": 0.240, "Cdl": 0.22, "sigma": 0.075},
        "unit_area": "Ω (normalised to electrode area)",
        "notes": "NCA // graphite cylindrical cell, 25°C, SOC=50%",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Polarization curve parameters
# ─────────────────────────────────────────────────────────────────────────────

POLARIZATION_PARAMS = {
    "PEM Fuel Cell — Pt/C (fresh)": {
        "type": "pem",
        "T":     353.15,   # K (80°C)
        "P_H2":  1.0,      # atm
        "P_O2":  0.21,     # atm (air)
        "i0":    8e-5,     # A cm⁻² (ORR exchange current density)
        "alpha": 0.50,     # transfer coefficient
        "R_mem": 0.16,     # Ω cm² (Nafion 212 hydrated)
        "i_L":   1.60,     # A cm⁻² (limiting current density)
        "B":     0.048,    # V (concentration loss coefficient)
    },
    "PEM Fuel Cell — Pt/C (degraded)": {
        "type": "pem",
        "T":     353.15,
        "P_H2":  1.0,
        "P_O2":  0.21,
        "i0":    2e-5,     # ECSA loss: ~75% reduction
        "alpha": 0.46,
        "R_mem": 0.38,     # membrane thinning + ionomer degradation
        "i_L":   1.10,     # flooding / GDL degradation
        "B":     0.058,
    },
    "NMC811 // Graphite": {
        "type": "battery",
        "cathode": "NMC811",
        "Q_cell":  195.0,  # mAh/g
        "C_rate":  1.0,
        "R_int":   0.018,  # Ω
        "i0":      4.5,    # mA/g
        "alpha":   0.5,
        "T":       298.15,
    },
    "LFP // Graphite": {
        "type": "battery",
        "cathode": "LFP",
        "Q_cell":  160.0,  # mAh/g
        "C_rate":  1.0,
        "R_int":   0.012,  # Ω
        "i0":      6.0,    # mA/g
        "alpha":   0.5,
        "T":       298.15,
    },
    "NCA // Graphite": {
        "type": "battery",
        "cathode": "NMC811",  # reuse NMC model; offset V handled in label
        "Q_cell":  205.0,
        "C_rate":  1.0,
        "R_int":   0.016,
        "i0":      5.0,
        "alpha":   0.5,
        "T":       298.15,
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Capacity fade parameters
# ─────────────────────────────────────────────────────────────────────────────

CAPACITY_FADE_PARAMS = {
    "NMC811 // Graphite": {
        "Q0":          195.0,    # mAh/g (initial capacity)
        "k_cal":       2.5e-3,   # calendar fade rate [day^-beta]
        "k_cyc":       6.0e-5,   # power-law cycle fade
        "beta_cal":    0.50,     # calendar exponent (SEI-diffusion limited)
        "beta_cyc":    0.55,     # cycle exponent
        "k_knee":      8.0e-5,   # Stage-2 fade coefficient
        "n_knee":      800,      # knee-point cycle number
        "Ea_cal":      24500,    # J/mol (calendar Arrhenius activation energy)
        "Ea_cyc":      28000,    # J/mol (cycle Arrhenius activation energy)
        "B_wang":      38.0,     # Wang model pre-factor
        "gamma":       1.25,     # C-rate exponent
        "R0":          0.018,    # Initial DC resistance [Ω]
        "k_sei":       4.0e-4,   # SEI resistance growth rate
        "end_of_life": 0.80,     # EOL capacity retention
    },
    "LFP // Graphite": {
        "Q0":          160.0,
        "k_cal":       1.2e-3,
        "k_cyc":       3.5e-5,
        "beta_cal":    0.50,
        "beta_cyc":    0.58,
        "k_knee":      4.0e-5,
        "n_knee":      1800,
        "Ea_cal":      22000,
        "Ea_cyc":      25000,
        "B_wang":      22.0,
        "gamma":       1.15,
        "R0":          0.012,
        "k_sei":       2.0e-4,
        "end_of_life": 0.80,
    },
    "NCA // Graphite": {
        "Q0":          205.0,
        "k_cal":       3.0e-3,
        "k_cyc":       7.5e-5,
        "beta_cal":    0.50,
        "beta_cyc":    0.52,
        "k_knee":      1.2e-4,
        "n_knee":      600,
        "Ea_cal":      27000,
        "Ea_cyc":      31000,
        "B_wang":      45.0,
        "gamma":       1.30,
        "R0":          0.016,
        "k_sei":       5.0e-4,
        "end_of_life": 0.80,
    },
    "PEM Fuel Cell — Pt/C": {
        "Q0":          1.0,          # normalised beginning-of-life power
        "k_deg":       6.5,          # μV/h steady-state voltage decay
        "k_startup":   3.2,          # μV per startup (load cycling)
        "ECSA0":       60.0,         # m²/g_Pt fresh
        "k_diss":      4.0e-5,       # Pt dissolution rate [m²/g per cycle]
        "k_osten":     2.5e-5,       # Ostwald ripening coefficient
        "R0_mem":      0.16,         # Ω·cm² fresh membrane resistance
        "k_mem":       2.0e-6,       # membrane drift [Ω·cm²/h]
        "k_corr":      1.0e-4,       # corrosion [Ω·cm²/h^0.5]
        "end_of_life": 0.90,         # EOL at 10% voltage loss
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Material summary for comparison radar chart
# Metrics: specific_capacity, rate_capability, calendar_stability,
#          cycle_stability, cost_index, energy_density
# All normalised to [0,1] scale relative to best-in-class.
# ─────────────────────────────────────────────────────────────────────────────

MATERIAL_COMPARISON = {
    "NMC811": {
        "Specific Capacity": 0.92,
        "Rate Capability":   0.75,
        "Calendar Stability": 0.58,
        "Cycle Stability":   0.62,
        "Cost Index":        0.55,   # lower is cheaper → inverted for display
        "Energy Density":    0.90,
        "Thermal Safety":    0.45,
    },
    "NMC622": {
        "Specific Capacity": 0.82,
        "Rate Capability":   0.78,
        "Calendar Stability": 0.68,
        "Cycle Stability":   0.72,
        "Cost Index":        0.65,
        "Energy Density":    0.80,
        "Thermal Safety":    0.58,
    },
    "LFP": {
        "Specific Capacity": 0.55,
        "Rate Capability":   0.85,
        "Calendar Stability": 0.90,
        "Cycle Stability":   0.95,
        "Cost Index":        0.90,
        "Energy Density":    0.52,
        "Thermal Safety":    0.95,
    },
    "NCA": {
        "Specific Capacity": 0.98,
        "Rate Capability":   0.72,
        "Calendar Stability": 0.52,
        "Cycle Stability":   0.55,
        "Cost Index":        0.48,
        "Energy Density":    0.95,
        "Thermal Safety":    0.38,
    },
    "Graphite Anode": {
        "Specific Capacity": 0.62,
        "Rate Capability":   0.65,
        "Calendar Stability": 0.80,
        "Cycle Stability":   0.85,
        "Cost Index":        0.95,
        "Energy Density":    0.60,
        "Thermal Safety":    0.78,
    },
    "Silicon Anode": {
        "Specific Capacity": 0.99,
        "Rate Capability":   0.45,
        "Calendar Stability": 0.42,
        "Cycle Stability":   0.35,
        "Cost Index":        0.50,
        "Energy Density":    0.98,
        "Thermal Safety":    0.65,
    },
    "Pt/C (PEM)": {
        "Specific Capacity": 0.85,
        "Rate Capability":   0.90,
        "Calendar Stability": 0.62,
        "Cycle Stability":   0.58,
        "Cost Index":        0.20,
        "Energy Density":    0.78,
        "Thermal Safety":    0.85,
    },
}
