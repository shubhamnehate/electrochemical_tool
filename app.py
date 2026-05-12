# Electrochemical Degradation & Performance Modelling Tool
# Run: streamlit run app.py

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

# Project modules
from models.eis_model import (
    generate_frequencies, generate_eis_data, randles_impedance,
    two_rc_warburg_impedance, pem_impedance, interpolate_params,
    fit_randles, fit_two_rc_warburg, compute_drt_approx,
)
from models.polarization import (
    pem_polarization_curve, pem_degradation_params,
    battery_voltage, battery_degradation_params, compute_dqdv,
)
from models.capacity_fade import (
    calendar_fade_sqrt, cycle_fade_power, cycle_fade_two_stage,
    combined_fade_wang, sei_resistance_growth,
    pem_voltage_decay, pem_ecsa_loss,
    effective_fade_rate, temperature_acceleration_factor,
)
from data.material_params import (
    EIS_PARAMS, POLARIZATION_PARAMS, CAPACITY_FADE_PARAMS, MATERIAL_COMPARISON,
)

# ── Page config ────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Electrochemical Degradation & Performance Modelling",
    page_icon="⚗",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Colour palette ─────────────────────────────────────────────────────────────
COLORS = {
    "fresh":      "#2ecc71",
    "aged":       "#e74c3c",
    "mid":        "#f39c12",
    "fit":        "#3498db",
    "secondary":  "#9b59b6",
    "bg":         "#0e1117",
    "grid":       "#2c2c2c",
}
MATERIAL_COLORS = px.colors.qualitative.Plotly

PLOT_LAYOUT = dict(
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(20,20,30,0.9)",
    font=dict(color="#e0e0e0", family="Inter, sans-serif"),
    xaxis=dict(gridcolor="#2a2a3a", zerolinecolor="#444"),
    yaxis=dict(gridcolor="#2a2a3a", zerolinecolor="#444"),
    legend=dict(bgcolor="rgba(30,30,40,0.8)", bordercolor="#444", borderwidth=1),
    margin=dict(l=60, r=30, t=50, b=60),
)

# ── Sidebar ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("⚗ Electrochemical Tool")
    st.caption("Degradation & Performance Modelling")
    st.divider()
    st.markdown("""
    **Navigation**
    Use the tabs to explore:
    - **EIS** — Nyquist / Bode / DRT / circuit fitting
    - **Polarization** — V-I curves, overpotentials, power
    - **Capacity Fade** — ageing models, resistance growth
    - **Comparison** — multi-material radar & bar charts
    """)
    st.divider()
    st.markdown("**Global settings**")
    noise_level = st.slider("Measurement noise (%)", 0, 10, 3, 1) / 100.0
    show_fit = st.checkbox("Show equivalent-circuit fit", value=True)
    dark_theme = True  # always dark for presentation

# ── Main tabs ──────────────────────────────────────────────────────────────────
tab_eis, tab_polar, tab_fade, tab_compare = st.tabs([
    "🔬 EIS Spectroscopy",
    "⚡ Polarization Curves",
    "📉 Capacity Fade",
    "📊 Material Comparison",
])


# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 – EIS SPECTROSCOPY
# ══════════════════════════════════════════════════════════════════════════════
with tab_eis:
    st.subheader("Electrochemical Impedance Spectroscopy")
    st.caption("Nyquist plots, Bode diagrams, DRT analysis, and equivalent circuit fitting")

    col_ctrl, col_plot = st.columns([1, 3], gap="large")

    with col_ctrl:
        eis_material = st.selectbox("Electrode material", list(EIS_PARAMS.keys()), key="eis_mat")
        deg_level = st.slider("Degradation level", 0.0, 1.0, 0.3, 0.01,
                               format="%.2f", key="eis_deg",
                               help="0 = fresh, 1 = fully degraded (literature end-of-life)")

        freq_low  = st.select_slider("Low frequency [Hz]",
                                      options=[0.001, 0.01, 0.1], value=0.01, key="f_low")
        freq_high = st.select_slider("High frequency [Hz]",
                                      options=[1e3, 1e4, 1e5], value=1e5, key="f_high")
        n_pts = st.slider("Number of frequency points", 20, 80, 50, 5, key="n_pts")

        view_mode = st.radio("Plot view", ["Nyquist", "Bode", "DRT", "All"],
                              index=0, key="eis_view")

        st.divider()
        st.markdown("**Interpolated circuit parameters**")

        mat = EIS_PARAMS[eis_material]
        circuit = mat["circuit"]
        params_interp = interpolate_params(mat["fresh"], mat["aged"], deg_level)

        for k, v in params_interp.items():
            st.metric(k, f"{v:.4f}", delta=None)

    with col_plot:
        freqs = generate_frequencies(freq_low, freq_high, n_pts)
        omega = 2 * np.pi * freqs

        # Generate synthetic EIS data
        if circuit == "randles":
            model_fn = lambda om, **p: randles_impedance(om, p["Rs"], p["Rct"], p["Cdl"], p["sigma"])
            Z_data = generate_eis_data(
                lambda om, **p: randles_impedance(om, **p),
                params_interp, freqs, noise_level=noise_level
            )
        elif circuit == "two_rc_warburg":
            def _two_rc(om, **p):
                return two_rc_warburg_impedance(om, p["Rs"], p["R_SEI"], p["C_SEI"],
                                                 p["Rct"], p["Cdl"], p["sigma"])
            Z_data = generate_eis_data(_two_rc, params_interp, freqs, noise_level=noise_level)
        else:  # pem
            def _pem(om, **p):
                return pem_impedance(om, p.get("Rmem", p.get("Rs", 0.14)),
                                      p["Rct"], p["Cdl"], p["sigma"], p["L"])
            Z_data = generate_eis_data(_pem, params_interp, freqs, noise_level=noise_level)

        # Fitting
        fit_result = None
        if show_fit and circuit in ("randles", "two_rc_warburg"):
            if circuit == "randles":
                fit_result = fit_randles(freqs, Z_data,
                                          p0=[params_interp["Rs"],
                                              params_interp["Rct"],
                                              params_interp["Cdl"],
                                              params_interp["sigma"]])
            else:
                p0 = [params_interp["Rs"], params_interp["R_SEI"], params_interp["C_SEI"],
                      params_interp["Rct"], params_interp["Cdl"], params_interp["sigma"]]
                fit_result = fit_two_rc_warburg(freqs, Z_data, p0=p0)

        # ── Plot builder ───────────────────────────────────────────────────────
        def make_nyquist():
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=Z_data.real, y=-Z_data.imag,
                mode="markers", name="Simulated EIS data",
                marker=dict(color=COLORS["fresh"], size=7, symbol="circle"),
            ))
            if fit_result and fit_result.get("success") and "Z_fit" in fit_result:
                Zf = fit_result["Z_fit"]
                fig.add_trace(go.Scatter(
                    x=Zf.real, y=-Zf.imag,
                    mode="lines", name=f"Circuit fit (err={fit_result['residual_pct']:.1f}%)",
                    line=dict(color=COLORS["fit"], width=2.5, dash="solid"),
                ))
            # Frequency labels at decades
            decade_idx = np.round(np.linspace(0, len(freqs)-1, 6)).astype(int)
            for i in decade_idx:
                fig.add_annotation(
                    x=Z_data[i].real, y=-Z_data[i].imag,
                    text=f"{freqs[i]:.2g} Hz",
                    showarrow=False, yshift=12,
                    font=dict(size=9, color="#aaaaaa"),
                )
            fig.update_layout(
                **PLOT_LAYOUT,
                title=f"Nyquist Plot — {eis_material} (deg={deg_level:.0%})",
                xaxis_title="Z' [Ω]",
                yaxis_title="-Z'' [Ω]",
                yaxis_scaleanchor="x",
            )
            return fig

        def make_bode():
            fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                                subplot_titles=["Magnitude |Z| [Ω]", "Phase angle [°]"])
            Z_mag   = np.abs(Z_data)
            Z_phase = np.degrees(np.angle(Z_data))

            fig.add_trace(go.Scatter(x=freqs, y=Z_mag, mode="markers+lines",
                                      name="|Z|", line=dict(color=COLORS["fresh"])),
                          row=1, col=1)
            fig.add_trace(go.Scatter(x=freqs, y=Z_phase, mode="markers+lines",
                                      name="Phase", line=dict(color=COLORS["secondary"])),
                          row=2, col=1)
            if fit_result and fit_result.get("success") and "Z_fit" in fit_result:
                Zf = fit_result["Z_fit"]
                fig.add_trace(go.Scatter(x=freqs, y=np.abs(Zf), mode="lines",
                                          name="Fit |Z|",
                                          line=dict(color=COLORS["fit"], dash="dash")),
                              row=1, col=1)
                fig.add_trace(go.Scatter(x=freqs, y=np.degrees(np.angle(Zf)), mode="lines",
                                          name="Fit phase",
                                          line=dict(color=COLORS["fit"], dash="dash")),
                              row=2, col=1)
            fig.update_xaxes(type="log", title_text="Frequency [Hz]", row=2, col=1,
                              gridcolor="#2a2a3a")
            fig.update_xaxes(type="log", gridcolor="#2a2a3a", row=1, col=1)
            fig.update_yaxes(gridcolor="#2a2a3a")
            fig.update_layout(**PLOT_LAYOUT,
                               title=f"Bode Plot — {eis_material} (deg={deg_level:.0%})")
            return fig

        def make_drt():
            log_tau, drt = compute_drt_approx(freqs, Z_data)
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=log_tau, y=drt, mode="lines", fill="tozeroy",
                name="DRT",
                line=dict(color=COLORS["mid"], width=2),
                fillcolor="rgba(243,156,18,0.25)",
            ))
            fig.update_layout(
                **PLOT_LAYOUT,
                title=f"Distribution of Relaxation Times — {eis_material}",
                xaxis_title="log₁₀(τ) [s]",
                yaxis_title="γ(τ) [normalised]",
            )
            return fig

        if view_mode == "Nyquist":
            st.plotly_chart(make_nyquist(), use_container_width=True)
        elif view_mode == "Bode":
            st.plotly_chart(make_bode(), use_container_width=True)
        elif view_mode == "DRT":
            st.plotly_chart(make_drt(), use_container_width=True)
        else:  # All
            c1, c2 = st.columns(2)
            with c1:
                st.plotly_chart(make_nyquist(), use_container_width=True)
                st.plotly_chart(make_drt(), use_container_width=True)
            with c2:
                st.plotly_chart(make_bode(), use_container_width=True)

        # Fit summary table
        if show_fit and fit_result and fit_result.get("success"):
            st.subheader("Equivalent circuit fit results")
            param_keys = [k for k in fit_result if not k.endswith("_err")
                          and k not in ("success", "Z_fit", "residual_pct")]
            rows = []
            for pk in param_keys:
                err_key = pk + "_err"
                rows.append({
                    "Parameter": pk,
                    "Fitted value": f"{fit_result[pk]:.5f}",
                    "Std. error":   f"±{fit_result.get(err_key, 0):.5f}",
                    "Literature (interpolated)": f"{params_interp.get(pk, '—'):.5f}"
                    if pk in params_interp else "—",
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
            st.caption(f"Mean relative residual: **{fit_result['residual_pct']:.2f}%**   |   "
                       f"Circuit: **{circuit}**   |   "
                       f"Notes: {mat['notes']}")


# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 – POLARIZATION CURVES
# ══════════════════════════════════════════════════════════════════════════════
with tab_polar:
    st.subheader("Polarization Curves")
    st.caption("V-I curves, power density, and individual overpotential contributions")

    col_c, col_p = st.columns([1, 3], gap="large")

    with col_c:
        cell_type = st.selectbox(
            "Cell / material", list(POLARIZATION_PARAMS.keys()), key="pol_mat"
        )
        pol_deg = st.slider("Degradation level", 0.0, 1.0, 0.0, 0.05, key="pol_deg",
                             format="%.2f")
        pol_temp = st.slider("Temperature [°C]", 15, 90, 25 if "Battery" not in cell_type else 25,
                              1, key="pol_temp")
        if POLARIZATION_PARAMS[cell_type]["type"] == "battery":
            c_rate = st.select_slider("C-rate", options=[0.1, 0.2, 0.5, 1.0, 2.0, 5.0],
                                       value=1.0, key="pol_crate")
        else:
            rh = st.slider("Relative humidity [%]", 50, 100, 100, 5, key="pol_rh")
            c_rate = 1.0

        show_contributions = st.checkbox("Show overpotential breakdown", value=True)
        show_power = st.checkbox("Show power density", value=True)

    with col_p:
        p = dict(POLARIZATION_PARAMS[cell_type])
        p["T"] = pol_temp + 273.15

        ptype = p["type"]

        if ptype == "pem":
            p_fresh = dict(p)
            p_deg = pem_degradation_params(p, pol_deg)
            i_array = np.linspace(0.001, p_deg["i_L"] * 0.97, 300)

            res_fresh = pem_polarization_curve(i_array, p_fresh)
            res_deg   = pem_polarization_curve(i_array, p_deg)

            nrows = 2 if show_power else 1
            fig = make_subplots(rows=nrows, cols=1, shared_xaxes=True, vertical_spacing=0.10,
                                subplot_titles=(["V-I Curve", "Power Density"] if show_power
                                                else ["V-I Curve"]))

            fig.add_trace(go.Scatter(x=i_array, y=res_fresh["V_cell"], name="Fresh",
                                      line=dict(color=COLORS["fresh"], width=2.5)),
                          row=1, col=1)
            fig.add_trace(go.Scatter(x=i_array, y=res_deg["V_cell"],
                                      name=f"Degraded ({pol_deg:.0%})",
                                      line=dict(color=COLORS["aged"], width=2.5)),
                          row=1, col=1)
            fig.add_hline(y=res_fresh["E_rev"][0], line=dict(color="#555", dash="dot"),
                          annotation_text=f"E_rev={res_fresh['E_rev'][0]:.3f} V",
                          annotation_position="right", row=1, col=1)

            if show_power:
                fig.add_trace(go.Scatter(x=i_array, y=res_fresh["P_density"] * 1000,
                                          name="P (fresh)", fill="tozeroy",
                                          line=dict(color=COLORS["fresh"]),
                                          fillcolor="rgba(46,204,113,0.15)"),
                              row=2, col=1)
                fig.add_trace(go.Scatter(x=i_array, y=res_deg["P_density"] * 1000,
                                          name=f"P (deg {pol_deg:.0%})",
                                          line=dict(color=COLORS["aged"]),
                                          fillcolor="rgba(231,76,60,0.15)", fill="tozeroy"),
                              row=2, col=1)
                fig.update_yaxes(title_text="Power density [mW cm⁻²]", row=2, col=1,
                                  gridcolor="#2a2a3a")
                fig.update_xaxes(title_text="Current density [A cm⁻²]", row=2, col=1,
                                  gridcolor="#2a2a3a")

            fig.update_yaxes(title_text="Cell voltage [V]", row=1, col=1, gridcolor="#2a2a3a")
            fig.update_xaxes(gridcolor="#2a2a3a", row=1, col=1)
            fig.update_layout(**PLOT_LAYOUT,
                               title=f"PEM Fuel Cell — {cell_type}")
            st.plotly_chart(fig, use_container_width=True)

            # Overpotential breakdown
            if show_contributions:
                st.subheader("Overpotential breakdown at selected operating point")
                op_i = st.slider("Operating current density [A cm⁻²]",
                                  float(i_array[1]), float(i_array[-5]),
                                  float(np.median(i_array)), key="op_i")
                idx = np.argmin(np.abs(i_array - op_i))
                labels = ["Activation η_act", "Ohmic η_ohm", "Concentration η_conc"]
                vals_f = [res_fresh["eta_act"][idx], res_fresh["eta_ohm"][idx],
                          res_fresh["eta_conc"][idx]]
                vals_d = [res_deg["eta_act"][idx], res_deg["eta_ohm"][idx],
                          res_deg["eta_conc"][idx]]
                bar_fig = go.Figure(data=[
                    go.Bar(name="Fresh", x=labels, y=vals_f,
                           marker_color=COLORS["fresh"]),
                    go.Bar(name=f"Degraded ({pol_deg:.0%})", x=labels, y=vals_d,
                           marker_color=COLORS["aged"]),
                ])
                bar_fig.update_layout(**PLOT_LAYOUT,
                                       title=f"Overpotentials at i={op_i:.3f} A cm⁻²",
                                       yaxis_title="Overpotential [V]",
                                       barmode="group")
                st.plotly_chart(bar_fig, use_container_width=True)

        else:  # battery
            p["C_rate"] = c_rate
            p_deg = battery_degradation_params(p, pol_deg)

            soc_array = np.linspace(0.98, 0.02, 300)
            res_fresh = battery_voltage(soc_array, p)
            res_deg   = battery_voltage(soc_array, p_deg)

            # Convert SOC to capacity
            Q_fresh = soc_array * p["Q_cell"]
            Q_deg   = soc_array * p_deg["Q_cell"]

            fig = make_subplots(rows=2 if show_power else 1, cols=1,
                                shared_xaxes=True, vertical_spacing=0.10,
                                subplot_titles=["Discharge Curve", "Power"] if show_power
                                               else ["Discharge Curve"])

            fig.add_trace(go.Scatter(x=Q_fresh, y=res_fresh["V_cell"], name="Fresh",
                                      line=dict(color=COLORS["fresh"], width=2.5)),
                          row=1, col=1)
            fig.add_trace(go.Scatter(x=Q_deg, y=res_deg["V_cell"],
                                      name=f"Degraded ({pol_deg:.0%})",
                                      line=dict(color=COLORS["aged"], width=2.5)),
                          row=1, col=1)
            fig.add_trace(go.Scatter(x=Q_fresh, y=res_fresh["V_ocv"], name="OCV",
                                      line=dict(color="#888", dash="dot", width=1.5)),
                          row=1, col=1)

            if show_power:
                fig.add_trace(go.Scatter(x=Q_fresh, y=res_fresh["P_cell"],
                                          name="P (fresh)", fill="tozeroy",
                                          line=dict(color=COLORS["fresh"]),
                                          fillcolor="rgba(46,204,113,0.15)"),
                              row=2, col=1)
                fig.add_trace(go.Scatter(x=Q_deg, y=res_deg["P_cell"],
                                          name=f"P (deg {pol_deg:.0%})",
                                          line=dict(color=COLORS["aged"]),
                                          fillcolor="rgba(231,76,60,0.15)", fill="tozeroy"),
                              row=2, col=1)
                fig.update_yaxes(title_text="Power [W g⁻¹]", row=2, col=1,
                                  gridcolor="#2a2a3a")
                fig.update_xaxes(title_text="Capacity [mAh g⁻¹]", row=2, col=1,
                                  gridcolor="#2a2a3a")

            fig.update_yaxes(title_text="Voltage [V]", row=1, col=1, gridcolor="#2a2a3a")
            fig.update_xaxes(gridcolor="#2a2a3a", row=1, col=1)
            fig.update_layout(**PLOT_LAYOUT,
                               title=f"Discharge Curve — {cell_type} @ {c_rate}C")
            st.plotly_chart(fig, use_container_width=True)

            # dQ/dV
            if show_contributions:
                st.subheader("dQ/dV Analysis")
                dcols = st.columns(2)
                for dcol, (res, label, col) in zip(dcols, [
                    (res_fresh, "Fresh", COLORS["fresh"]),
                    (res_deg, f"Degraded ({pol_deg:.0%})", COLORS["aged"]),
                ]):
                    soc_use = soc_array if label == "Fresh" else soc_array
                    q_use = (soc_use * (p["Q_cell"] if label == "Fresh" else p_deg["Q_cell"]))
                    V_mid, dqdv = compute_dqdv(
                        soc_use, res["V_cell"], p["Q_cell"] if label == "Fresh" else p_deg["Q_cell"]
                    )
                    dqdv_fig = go.Figure(go.Scatter(
                        x=V_mid, y=dqdv, mode="lines", fill="tozeroy",
                        line=dict(color=col, width=2),
                        fillcolor=col.replace(")", ",0.2)").replace("rgb", "rgba"),
                    ))
                    dqdv_fig.update_layout(
                        **PLOT_LAYOUT,
                        title=f"dQ/dV — {label}",
                        xaxis_title="Voltage [V]",
                        yaxis_title="dQ/dV [mAh g⁻¹ V⁻¹]",
                        height=280,
                    )
                    with dcol:
                        st.plotly_chart(dqdv_fig, use_container_width=True)


# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 – CAPACITY FADE
# ══════════════════════════════════════════════════════════════════════════════
with tab_fade:
    st.subheader("Capacity Fade & Lifetime Modelling")
    st.caption("Calendar ageing, cycle ageing, SEI growth, temperature and C-rate effects")

    col_cf, col_fp = st.columns([1, 3], gap="large")

    with col_cf:
        fade_mat = st.selectbox("Material / cell", list(CAPACITY_FADE_PARAMS.keys()), key="fd_mat")
        fp = CAPACITY_FADE_PARAMS[fade_mat]
        is_pem = fade_mat == "PEM Fuel Cell — Pt/C"

        if not is_pem:
            max_cycles = st.slider("Max cycles", 200, 3000, 1500, 100, key="fd_cyc")
            T_degC     = st.slider("Temperature [°C]", 15, 60, 25, 5, key="fd_T")
            crate_fade = st.select_slider("C-rate", [0.5, 1.0, 2.0, 3.0, 5.0],
                                           value=1.0, key="fd_crate")
            show_calendar = st.checkbox("Show calendar ageing", value=True)
            show_resistance = st.checkbox("Show resistance growth", value=True)
            fade_model = st.radio("Cycle fade model",
                                   ["Power-law", "Two-stage (knee)", "Wang et al. (combined)"],
                                   key="fd_model")
        else:
            max_hours = st.slider("Operating hours", 500, 20000, 8000, 500, key="fd_hrs")
            n_startups = st.slider("Load cycles (startups)", 0, 5000, 1000, 100, key="fd_st")

    with col_fp:
        if not is_pem:
            T_K = T_degC + 273.15
            n_cycles = np.arange(0, max_cycles + 1, 1, dtype=float)
            Q0 = fp["Q0"]

            # Effective fade rate with T and C-rate acceleration
            k_cyc_eff = effective_fade_rate(
                fp["k_cyc"], T_K, crate_fade,
                Ea=fp["Ea_cyc"], gamma=fp.get("gamma", 1.25)
            )

            if fade_model == "Power-law":
                Q_cyc = cycle_fade_power(n_cycles, Q0, k_cyc_eff, fp["beta_cyc"])
            elif fade_model == "Two-stage (knee)":
                Q_cyc = cycle_fade_two_stage(n_cycles, Q0, k_cyc_eff,
                                              fp["k_knee"], fp["n_knee"])
            else:
                Q_cyc = combined_fade_wang(n_cycles, n_cycles / 365.0, Q0,
                                            fp["B_wang"], fp["Ea_cyc"], T_K, crate_fade)

            # EOL marker
            eol_cap = fp["Q0"] * fp["end_of_life"]
            eol_idx = np.searchsorted(-Q_cyc, -eol_cap)
            eol_cycle = n_cycles[eol_idx] if eol_idx < len(n_cycles) else None

            # Calendar ageing (vs days)
            t_days = np.linspace(0, max_cycles * 2, 500)
            k_cal_eff = effective_fade_rate(fp["k_cal"], T_K, 1.0, Ea=fp["Ea_cal"])
            Q_cal = calendar_fade_sqrt(t_days, Q0, k_cal_eff, fp["beta_cal"])

            # Resistance growth
            R_growth = sei_resistance_growth(n_cycles, fp["R0"], fp["k_sei"])

            # ── Capacity retention vs cycles ──────────────────────────────────
            fig_cap = go.Figure()
            fig_cap.add_trace(go.Scatter(
                x=n_cycles, y=Q_cyc / Q0 * 100,
                mode="lines", name=fade_model,
                line=dict(color=COLORS["fresh"], width=2.5),
            ))
            # Reference: 1C, 25°C baseline
            Q_ref = cycle_fade_power(n_cycles, Q0, fp["k_cyc"], fp["beta_cyc"])
            fig_cap.add_trace(go.Scatter(
                x=n_cycles, y=Q_ref / Q0 * 100,
                mode="lines", name="Baseline (1C, 25°C)",
                line=dict(color="#666", width=1.5, dash="dot"),
            ))
            fig_cap.add_hline(y=80, line=dict(color=COLORS["aged"], dash="dash"),
                               annotation_text="EOL threshold (80%)",
                               annotation_font_color=COLORS["aged"])
            if eol_cycle:
                fig_cap.add_vline(x=eol_cycle, line=dict(color=COLORS["aged"], dash="dashdot"),
                                   annotation_text=f"EOL ~{int(eol_cycle)} cycles",
                                   annotation_font_color=COLORS["aged"])
            fig_cap.update_layout(
                **PLOT_LAYOUT,
                title=f"Capacity Retention — {fade_mat}   ({T_degC}°C, {crate_fade}C)",
                xaxis_title="Cycle number",
                yaxis_title="Capacity retention [%]",
                yaxis_range=[50, 102],
            )
            st.plotly_chart(fig_cap, use_container_width=True)

            # ── Multi-temperature comparison ──────────────────────────────────
            st.subheader("Temperature sensitivity")
            temp_compare_figs = go.Figure()
            temps_K = [273.15+15, 273.15+25, 273.15+35, 273.15+45, 273.15+60]
            temp_palette = ["#85c1e9", "#2ecc71", "#f1c40f", "#e67e22", "#e74c3c"]
            for T_i, col_i in zip(temps_K, temp_palette):
                k_i = effective_fade_rate(fp["k_cyc"], T_i, crate_fade,
                                           Ea=fp["Ea_cyc"], gamma=fp.get("gamma", 1.25))
                Q_i = cycle_fade_power(n_cycles, Q0, k_i, fp["beta_cyc"])
                temp_compare_figs.add_trace(go.Scatter(
                    x=n_cycles, y=Q_i / Q0 * 100,
                    mode="lines", name=f"{int(T_i-273.15)}°C",
                    line=dict(color=col_i, width=2),
                ))
            temp_compare_figs.add_hline(y=80, line=dict(color="#555", dash="dash"))
            temp_compare_figs.update_layout(
                **PLOT_LAYOUT,
                title="Capacity vs Cycle Number at Different Temperatures",
                xaxis_title="Cycle number",
                yaxis_title="Capacity retention [%]",
                yaxis_range=[50, 102],
            )
            st.plotly_chart(temp_compare_figs, use_container_width=True)

            if show_calendar:
                fig_cal = go.Figure()
                fig_cal.add_trace(go.Scatter(
                    x=t_days / 365, y=Q_cal / Q0 * 100,
                    mode="lines", name="Calendar fade (SEI growth)",
                    line=dict(color=COLORS["mid"], width=2.5), fill="tozeroy",
                    fillcolor="rgba(243,156,18,0.12)",
                ))
                fig_cal.add_hline(y=80, line=dict(color=COLORS["aged"], dash="dash"),
                                   annotation_text="EOL (80%)")
                fig_cal.update_layout(
                    **PLOT_LAYOUT,
                    title=f"Calendar Ageing — {fade_mat}   ({T_degC}°C, SOC=50%)",
                    xaxis_title="Storage time [years]",
                    yaxis_title="Capacity retention [%]",
                    yaxis_range=[60, 102],
                )
                st.plotly_chart(fig_cal, use_container_width=True)

            if show_resistance:
                fig_res = go.Figure()
                fig_res.add_trace(go.Scatter(
                    x=n_cycles, y=R_growth * 1000,
                    mode="lines", name="DC resistance (SEI growth)",
                    line=dict(color=COLORS["secondary"], width=2.5),
                ))
                fig_res.update_layout(
                    **PLOT_LAYOUT,
                    title=f"Resistance Growth — {fade_mat}",
                    xaxis_title="Cycle number",
                    yaxis_title="DC resistance [mΩ]",
                )
                st.plotly_chart(fig_res, use_container_width=True)

        else:
            # PEM fuel cell degradation
            hours = np.linspace(0, max_hours, 1000)
            V0 = 0.70   # V (BOL operating voltage at 0.6 A cm⁻²)

            V_deg = pem_voltage_decay(hours, V0, fp["k_deg"], fp["k_startup"],
                                       n_startups * hours / max_hours)
            n_cyc_pem = n_startups * hours / max_hours
            ECSA = pem_ecsa_loss(n_cyc_pem, fp["ECSA0"], fp["k_diss"], fp["k_osten"])

            from models.capacity_fade import ohmic_resistance_growth_pem
            R_pem = ohmic_resistance_growth_pem(hours, fp["R0_mem"], fp["k_mem"], fp["k_corr"])

            fig_pem = make_subplots(rows=3, cols=1, shared_xaxes=True,
                                     vertical_spacing=0.08,
                                     subplot_titles=["Cell voltage [V]",
                                                      "Pt ECSA [m² g⁻¹]",
                                                      "Membrane resistance [Ω cm²]"])
            fig_pem.add_trace(go.Scatter(x=hours/1000, y=V_deg,
                                          name="Cell voltage",
                                          line=dict(color=COLORS["fresh"], width=2)),
                              row=1, col=1)
            eol_V = V0 * fp["end_of_life"]
            fig_pem.add_hline(y=eol_V, line=dict(color=COLORS["aged"], dash="dash"),
                               annotation_text=f"EOL ({fp['end_of_life']:.0%} of BOL)",
                               row=1, col=1)
            fig_pem.add_trace(go.Scatter(x=hours/1000, y=ECSA,
                                          name="Pt ECSA",
                                          line=dict(color=COLORS["mid"], width=2)),
                              row=2, col=1)
            fig_pem.add_trace(go.Scatter(x=hours/1000, y=R_pem,
                                          name="R_mem",
                                          line=dict(color=COLORS["secondary"], width=2)),
                              row=3, col=1)
            fig_pem.update_xaxes(title_text="Operating time [kh]", row=3, col=1,
                                  gridcolor="#2a2a3a")
            for r in range(1, 4):
                fig_pem.update_xaxes(gridcolor="#2a2a3a", row=r, col=1)
                fig_pem.update_yaxes(gridcolor="#2a2a3a", row=r, col=1)
            fig_pem.update_layout(**PLOT_LAYOUT,
                                   title="PEM Fuel Cell Degradation Dashboard",
                                   height=650)
            st.plotly_chart(fig_pem, use_container_width=True)


# ══════════════════════════════════════════════════════════════════════════════
# TAB 4 – MATERIAL COMPARISON
# ══════════════════════════════════════════════════════════════════════════════
with tab_compare:
    st.subheader("Multi-Material Comparison")
    st.caption("Radar chart and bar chart comparison of key electrochemical performance metrics")

    col_m, col_chart = st.columns([1, 3], gap="large")

    with col_m:
        all_mats = list(MATERIAL_COMPARISON.keys())
        selected_mats = st.multiselect(
            "Select materials to compare",
            all_mats,
            default=["NMC811", "LFP", "NCA", "Graphite Anode"],
            key="cmp_mats",
        )
        chart_type = st.radio("Chart type", ["Radar", "Bar chart", "Both"], index=2,
                               key="cmp_chart")

        metrics = list(list(MATERIAL_COMPARISON.values())[0].keys())
        selected_metrics = st.multiselect("Metrics to display", metrics,
                                           default=metrics, key="cmp_metrics")

    with col_chart:
        if not selected_mats or not selected_metrics:
            st.info("Select at least one material and one metric to display.")
        else:
            def make_radar():
                fig = go.Figure()
                for idx, mat_name in enumerate(selected_mats):
                    vals = [MATERIAL_COMPARISON[mat_name].get(m, 0) for m in selected_metrics]
                    vals_closed = vals + [vals[0]]  # close the loop
                    fig.add_trace(go.Scatterpolar(
                        r=vals_closed,
                        theta=selected_metrics + [selected_metrics[0]],
                        fill="toself",
                        name=mat_name,
                        line=dict(color=MATERIAL_COLORS[idx % len(MATERIAL_COLORS)]),
                        fillcolor=MATERIAL_COLORS[idx % len(MATERIAL_COLORS)].replace(
                            "rgb", "rgba").replace(")", ",0.15)"),
                    ))
                fig.update_layout(
                    polar=dict(
                        radialaxis=dict(visible=True, range=[0, 1.05],
                                        tickfont=dict(size=9, color="#aaa"),
                                        gridcolor="#333"),
                        angularaxis=dict(tickfont=dict(size=11, color="#ddd"),
                                          gridcolor="#333"),
                        bgcolor="rgba(20,20,30,0.9)",
                    ),
                    paper_bgcolor="rgba(0,0,0,0)",
                    legend=dict(bgcolor="rgba(30,30,40,0.8)", bordercolor="#444",
                                borderwidth=1, font=dict(color="#ddd")),
                    title=dict(text="Normalised Performance Metrics — Radar",
                               font=dict(color="#e0e0e0")),
                    margin=dict(l=60, r=60, t=80, b=60),
                )
                return fig

            def make_bar():
                rows = []
                for mat_name in selected_mats:
                    for metric in selected_metrics:
                        rows.append({
                            "Material": mat_name,
                            "Metric": metric,
                            "Score": MATERIAL_COMPARISON[mat_name].get(metric, 0),
                        })
                df = pd.DataFrame(rows)
                fig = px.bar(df, x="Metric", y="Score", color="Material",
                             barmode="group",
                             color_discrete_sequence=MATERIAL_COLORS,
                             template="plotly_dark")
                fig.update_layout(
                    **PLOT_LAYOUT,
                    title="Normalised Performance Metrics — Bar Chart",
                    yaxis_range=[0, 1.1],
                    xaxis_tickangle=-30,
                )
                return fig

            if chart_type == "Radar":
                st.plotly_chart(make_radar(), use_container_width=True)
            elif chart_type == "Bar chart":
                st.plotly_chart(make_bar(), use_container_width=True)
            else:
                c1, c2 = st.columns(2)
                with c1:
                    st.plotly_chart(make_radar(), use_container_width=True)
                with c2:
                    st.plotly_chart(make_bar(), use_container_width=True)

            # Numeric comparison table
            st.subheader("Numeric scores (normalised 0–1)")
            df_tbl = pd.DataFrame(
                {mat: {m: MATERIAL_COMPARISON[mat].get(m, 0) for m in selected_metrics}
                 for mat in selected_mats}
            ).T.round(3)
            st.dataframe(df_tbl.style.background_gradient(cmap="RdYlGn", axis=None),
                         use_container_width=True)

            st.caption("Scores normalised to best-in-class (1.0) across published literature values. "
                       "Cost Index: higher = lower cost. "
                       "Thermal Safety: higher = safer.")


st.divider()
