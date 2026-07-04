# Rotating Electro-Optic / Geometric-Phase Frequency Shifter Modeling Handoff

This package is a Jupyter-first modeling handoff for the bulk z-cut LiNbO3 rotating-electro-optic project.

The goal is to get to an "80% of Ansys" level of design confidence without using a full commercial EM/RF simulator. The planned model stack is:

1. **Electro-optic tensor model** for z-cut LiNbO3 propagation along the optic axis.
2. **2D/2.5D electrostatic model** of the four-electrode geometry to compute field-per-volt, field uniformity, and the capacitance matrix.
3. **RF load/current model** for the EO electrode capacitance, including high-current limits at 6-100+ MHz.
4. **LDMOS/GaN RF PA output-network model** where the EO capacitance is part of the final load/match rather than a random external load.
5. **Quadrature-error and optical-efficiency model** to connect amplitude/phase errors to rotating eigenaxis errors, spurs, and geometric-phase frequency-shift performance.
6. **Predistortion / closed-loop calibration scaffold**, eventually using measured S-parameters and electrode-voltage pickups to generate corrected drive waveforms.

## Suggested workflow

Start here:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter lab
```

Then work through notebooks in order:

```text
00_project_constants.ipynb
01_eo_retardance_and_vpi.ipynb
02_electrostatic_field_and_capacitance.ipynb
03_rf_load_current_budget.ipynb
04_ldmos_output_network_screening.ipynb
05_quadrature_chirp_and_dpd_scaffold.ipynb
```

## Most important design question

For one EO electrode pair, the current requirement is not optional:

```math
i_C(t)=C_e\frac{dV_{diff}}{dt}
```

For a sinusoid:

```math
I_{peak}=2\pi f C_e V_{diff,peak}
```

For the present rough full-face geometry, use `C_e ≈ 13 pF` until measured otherwise. At `100 MHz`, `480 Vpeak differential`, this is about `3.9 Apeak` at the EO electrodes. That number drives the RF architecture.

## Design stance captured in this handoff

For broadband chirps, avoid making the core architecture a high-Q resonator bank. Instead, prefer a broadband class-AB-ish RF PA architecture:

```text
AWG/DDS chirp source
  -> driver / 0°/180° phase splitter
  -> push-pull LDMOS or GaN RF PA
  -> custom output network including the EO capacitance
  -> EO electrode pair
  -> capacitive HV pickup
  -> calibration / predistortion
```

The output network is not a generic 50-ohm output plus a random external capacitance. The EO capacitance is part of the output design.

