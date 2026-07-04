# Jupyter Modeling Plan

## Philosophy

The target is not Ansys-level field fidelity. The target is a design-grade notebook suite that can tell us:

- what voltage is actually needed,
- what current the RF stage must source,
- how electrode geometry changes field-per-volt and capacitance,
- what output network families are plausible,
- and how quadrature errors translate to optical spurs.

## Notebook 00: constants and geometry

Define all project constants in one place:

- wavelength
- LiNbO3 indices and EO coefficients
- crystal dimensions
- intended cut and beam direction
- electrode gap and geometry assumptions
- capacitance assumption / measured values
- frequency range
- voltage target

Output: a design table used by all notebooks.

## Notebook 01: EO retardance and Vpi

Implement:

- transverse EO impermeability matrix
- eigenvalues/eigenvectors
- retardance
- half-wave field
- Vpi from field-per-volt
- sensitivity sweeps vs `L`, `d`, `r22`, `n_o`, wavelength

Output: optical/electro-optic spec sheet.

## Notebook 02: electrostatic field and capacitance

Start with 2D cross-section electrostatics.

Solve:

```math
\nabla\cdot(\epsilon\nabla V)=0
```

Required outputs:

- field map `Ex(x,y), Ey(x,y)`
- field at beam center
- Gaussian-beam-weighted average field
- field uniformity across beam
- capacitance matrix of electrodes
- x/y electrode cross-coupling

Implementation options:

- first pass: finite difference grid for robustness
- better pass: `scikit-fem` mesh with dielectric regions and electrode boundaries

## Notebook 03: RF load/current budget

Use capacitance matrix or scalar capacitance to compute:

- `Z(f)`
- `Ipeak(f)`
- reactive apparent power
- dielectric loss estimate using tan-delta
- impact of extra parasitic pF
- transformer reflection scaling
- amp current target vs turns ratio

Output: RF current specification.

## Notebook 04: LDMOS output-network screening

Build ABCD/scikit-rf network models for candidate output networks:

- ideal voltage source with source resistance
- series L/R
- shunt R
- low-Q peaking
- transformer model with leakage and parasitic capacitance
- EO capacitance as final load

Optimize network values using `scipy.optimize.least_squares` or `differential_evolution`.

Objective examples:

- maximize minimum `|V_EO|` over band
- minimize phase ripple after group delay removal
- constrain current and resistor power
- penalize sharp resonances / high Q

## Notebook 05: quadrature chirp and DPD scaffold

Given two measured/modelled transfer functions:

```math
H_x(f), H_y(f)
```

compute drive predistortion:

```math
V_{drive,x}(f) = V_{target,x}(f) / H_x(f)
```

Then simulate:

- actual `Ex(t), Ey(t)`
- eigenaxis angle `theta(t)`
- retardance `Gamma(t)`
- Jones matrix evolution
- output optical sidebands/spurs

Start with linear inverse filtering. Later extend to ML-DPD using measured feedback.

