# Prompt to Give Claude

You are helping build a Jupyter/Python modeling package for a rotating electro-optic / geometric-phase optical frequency shifter based on a bulk z-cut LiNbO3 crystal.

Read the files in this package first, especially:

- `README.md`
- `docs/project_context.md`
- `docs/eo_model.md`
- `docs/rf_design_problem.md`
- `docs/modeling_plan.md`
- `docs/source_map.md`

Then implement the notebook suite and Python package described below.

## Core goals

We want an "80% of Ansys" Jupyter workflow, not just hand calculations. The package should calculate:

1. z-cut LiNbO3 EO retardance and Vpi using the impermeability tensor and eigenvalue splitting.
2. Electrostatic field-per-volt and capacitance matrix of a four-electrode cross-section.
3. RF current and apparent-power burden for a 13 pF-ish EO electrode load from 6-100+ MHz.
4. Candidate broadband LDMOS/GaN output-network behavior where the EO capacitance is part of the final output network.
5. Quadrature-drive amplitude/phase errors, optical eigenaxis rotation, and geometric-phase waveform errors.
6. Predistortion / closed-loop calibration scaffold using measured transfer functions or S-parameters.

## Hard constraints / design stance

- Do not model the RF stage as a generic 50-ohm amplifier followed by a random 13 pF load.
- The EO capacitance must be included in the final output-network model.
- Avoid high-Q resonator banks as the default because the desired application includes broadband chirps.
- The notebooks should make the high-current burden explicit through `i=C*dV/dt`.
- The model must separate electrode-voltage target, RF current, real power dissipation, and apparent reactive power.
- Use SI units internally.
- Every function should have docstrings and units in variable names or comments.
- Include plots and summary tables.
- Write tests for core formulae.

## Preferred packages

Use:

- numpy, scipy, matplotlib, pandas
- scikit-fem for FEM electrostatics
- scikit-rf for S-parameter/network modeling
- scikit-learn for simple DPD models later
- sympy where symbolic checks help

## Implementation order

1. Make `src/rot_eo_model/eo.py` robust and tested.
2. Make `src/rot_eo_model/rf_lumped.py` robust and tested.
3. Implement notebook 01 and 03 first because they are decisive.
4. Add a finite-difference electrostatics fallback before the scikit-fem version, so the notebook runs on any machine.
5. Add scikit-rf network model and optimizer.
6. Add quadrature/DPD scaffold.

## Do not overclaim

If a model is quasi-static, say so. If a model ignores transmission-line effects, say so. If it is not a substitute for HFSS/Ansys, say so. The point is to set design specs and kill bad architectures early.

