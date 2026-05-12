# Electrochemical Degradation & Performance Modelling Tool

**Personal Project — 2024**

An interactive simulation tool built with Python and Streamlit for modelling electrochemical degradation and performance across PEM fuel cell and Li-ion battery electrode materials. Parameters are derived from published degradation datasets and fitted against standard equivalent circuit models.

## Features

- **EIS Spectroscopy** — simulate and fit Nyquist plots, Bode diagrams, and Distribution of Relaxation Times (DRT) for Randles, two-RC-Warburg, and PEM equivalent circuits across selectable frequency ranges
- **Polarization Curves** — model V-I curves and power density for PEM fuel cells and Li-ion cells (NMC811, LFP, NCA); decompose activation, ohmic, and concentration overpotentials; overlay fresh vs degraded states; visualise dQ/dV profiles
- **Capacity Fade Dynamics** — simulate calendar ageing (SEI-diffusion model), cycle ageing (power-law, two-stage knee-point, Wang combined model), Arrhenius temperature acceleration, C-rate sensitivity, and resistance growth over lifetime
- **Material Comparison** — radar chart and grouped bar chart comparison of up to 7 electrode materials across 7 normalised performance metrics for rapid selection and benchmarking

## Materials Covered

| Material | Type |
|----------|------|
| NMC811 Cathode | Li-ion battery |
| LFP Cathode | Li-ion battery |
| NCA Cathode | Li-ion battery |
| Graphite Anode | Li-ion battery |
| Silicon Anode | Li-ion battery |
| Pt/C ORR Catalyst | PEM fuel cell |

## Tech Stack

- **Python 3.10+**
- **Streamlit** — interactive web interface
- **Plotly** — Nyquist, Bode, DRT, V-I, capacity fade, and radar visualisations
- **SciPy** — nonlinear least-squares equivalent circuit fitting (`curve_fit`)
- **NumPy / Pandas** — numerical computation and tabular display

## Project Structure

```
electrochemical_tool/
├── app.py                  # Main Streamlit application
├── models/
│   ├── eis_model.py        # EIS circuit impedance and fitting
│   ├── polarization.py     # PEM and battery V-I curve models
│   └── capacity_fade.py    # Ageing and degradation dynamics
├── data/
│   └── material_params.py  # Electrode material parameter database
└── requirements.txt
```

## Getting Started

```bash
pip install -r requirements.txt
streamlit run app.py
```

Open `http://localhost:8501` in your browser.

## Motivation

Built to support rapid material comparison and degradation analysis during electrode selection workflows — specifically to visualise how EIS spectra evolve with ageing, how polarisation performance degrades under different operating conditions, and how capacity retention varies across chemistries and temperatures without requiring experimental data for each scenario.
