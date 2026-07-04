# ADS / SPICE Import Package for Rotating EO Modulator

This directory contains the RF-ready linear load models for the locked 500 µm hybrid geometry.

> **WARNING**: These are strictly small-signal linear load models (capacitive with optional dielectric loss). They do **NOT** model the nonlinear LDMOS PA or the RF matching network itself.

## Touchstone Files (.sNp)
Use these in ADS using the `SnP` component.

### 4-Port Single-Ended (`eo_load_4port_single_ended_*.s4p`)
Reference impedance: 50 ohms.
Port definitions:
1. `+x_full` electrode (left face)
2. `-x_full` electrode (right face)
3. `+y_strip` electrode (top face)
4. `-y_strip` electrode (bottom face)

### 2-Port Differential (`eo_load_2port_diff_*.s2p`)
Reference impedance: 50 ohms or 100 ohms (as marked in filename).
Port definitions:
1. `x` differential pair (driven as $+V/2$ and $-V/2$)
2. `y` differential pair (driven as $+V/2$ and $-V/2$)

Frequency Range: 0.5 MHz to 200 MHz

## SPICE Subcircuits (.sp)
Standard text subcircuits defining pairwise Maxwell capacitances.
`EO_LOAD_4NODE p_x m_x p_y m_y`

For the lossy variants (`tand_1e-4` etc.), parallel resistors are added to simulate the equivalent dielectric conductance evaluated at **100 MHz**. For accurate broadband loss simulation in SPICE, consider using Laplace blocks or replacing the `.sp` with the corresponding Touchstone file in your SPICE simulator.
