# ADS / SPICE Import Package for Rotating EO Modulator

This directory contains the RF-ready linear load models for the locked 500 µm hybrid geometry.

> **WARNING**: These are strictly small-signal linear load models. They do **NOT** model the nonlinear LDMOS PA or the matching network.

## Touchstone Files (.sNp)
Use these in ADS using the `SnP` component. The files are written in **RI** (Real/Imaginary) format.
Includes `_tol_0.8` and `_tol_1.2` tolerance variants for the 4-port core.

## Minimal ADS SnP Smoke Test
To verify format compliance in ADS before proceeding to layout/matching:
1. Place an `SnP` component in an empty ADS schematic.
2. Point it to `eo_load_4port_single_ended_R50.s4p`.
3. Terminate all four ports with 50 Ohm (`Term` components).
4. Run an S-parameter sweep from 0.5 MHz to 200 MHz.
5. Check expected input capacitance (Y-parameters) at low frequency:
   - $C_{diff,x} \approx 16.88$ pF
   - $C_{diff,y} \approx 16.18$ pF
6. Check input impedance magnitude $|Z|$ at 100 MHz:
   - $|Z_{diff,x}| \approx 94\ \Omega$
   - $|Z_{diff,y}| \approx 98\ \Omega$

**IMPORTANT**: The 4-port `.s4p` is the definitive authoritative physical model. All other `.s2p` models are derivative reductions.

## SPICE Subcircuits
`eo_load_4node.sp` is the authoritative pure-capacitance Maxwell model.

**Important Note on Lossy Models**:
- The `_at_100MHz.sp` files contain fixed resistors valid EXACTLY at 100 MHz.
- The `_ladder_0p5MHz_200MHz.sub` files are broadband passive RC-ladder approximations that maintain a nearly constant $\tan\delta$ over the specified bandwidth. These are robust for transient and AC analysis.

## Fixture Parasitics Wrapper
The `fixture_parasitics_wrapper.net` (ADS syntax) demonstrates how to properly wrap the ideal 4-port core with estimated series inductances and resistances. *These parasitics must be verified against VNA measurements.*

## Measurement Fitting
`fit_vna_measurements.py` provides a scaffold to extract $C_{node}$, $\tan\delta$, and series $R/L$ directly from VNA `.s4p` files.
