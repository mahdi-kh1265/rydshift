# Project Context

## Project name / framing

Working project: **rotating electro-optic / geometric-phase optical waveform generator**.

Core concept: use bulk z-cut lithium niobate (LiNbO3) with RF quadrature transverse electric fields to create a time-dependent birefringent ellipse in the transverse plane. If the retardance is near half-wave and the eigenaxes rotate in time, circularly polarized input light can acquire a geometric/Pancharatnam-Berry phase whose time derivative produces an optical frequency shift or chirp.

Preferred high-level framing:

> RF-programmable, beam-stationary geometric-phase optical detuning/waveform generator for chirped Rydberg excitation.

The immediate modeling goal is not a full Ansys/HFSS/Maxwell simulation. The goal is a Jupyter/Python calculation suite that is much better than back-of-the-envelope and good enough to set design specs before hardware.

## Physical crystal assumption

Current working geometry:

- Material: lithium niobate, LiNbO3
- Cut: z-cut
- Dimensions under discussion: `3 mm x 3 mm x 30 mm`
- Intended optical path length: ideally `L = 30 mm` along crystal `z` / optic axis
- Transverse electrode gap scale: about `d = 3 mm`
- Wavelength: `lambda0 = 780 nm`
- Ordinary index: `n_o ~ 2.286` at 780 nm, to be refined with Sellmeier/vendor data
- EO coefficient used in rotating-EO derivation: `|r_22| ~ 6.8 pm/V`, to be refined by source and wavelength
- Rough full-face one-axis capacitance currently treated as `C_e ~ 13 pF`

**Critical geometry warning:** vendor language for `z-cut 3x3x30 mm` must be checked. If the z-axis is normal to a 3x30 face, then normal incidence may give a 3 mm optical length instead of 30 mm. The Vpi result changes by about 10x. The notebook should make this impossible to overlook.

## Why z-cut matters

For unperturbed LN in crystal axes:

```math
\boldsymbol{\eta}^{(0)} = \begin{pmatrix}
\eta_o & 0 & 0 \\
0 & \eta_o & 0 \\
0 & 0 & \eta_e
\end{pmatrix}
```

where:

```math
\eta_o = \frac{1}{n_o^2}, \quad \eta_e = \frac{1}{n_e^2}
```

For z-cut normal incidence, `k || z`, so the transverse optical plane is `x-y`:

```math
\boldsymbol{\eta}^{(0)}_\perp = \eta_o I
```

This means no static transverse birefringence: both transverse polarizations see `n_o`. The EO drive creates the birefringence and rotates the eigenaxes.

## Current preferred RF architecture

The present preferred architecture for broadband chirps is **not** a high-Q resonant tank as the core design. Resonators are excellent for single-frequency voltage enhancement, but fight chirps.

Preferred direction:

```text
AWG/DDS chirp source
  -> RF driver / phase splitter
  -> push-pull LDMOS or GaN RF PA, likely class-AB-ish
  -> output network designed around EO capacitance
  -> EO electrode pair
  -> capacitive HV pickup
  -> calibration / predistortion
```

Important architecture principle:

> The EO capacitance is part of the final RF output network. Do not design a generic 50-ohm amplifier and then bolt a 13 pF EO load onto it.

## Why the RF problem is hard

For a capacitive EO load:

```math
I_C = 2\pi f C_e V_{diff,peak}
```

With `C_e = 13 pF`, `Vdiff = 480 Vpeak`, `f = 100 MHz`:

```math
I_C \approx 3.9 Apeak
```

That is the current at the EO electrodes. Any transformer/output network can move this burden around, but cannot make it disappear.

## Working RF part candidates

Likely LDMOS starting candidates:

- NXP MRF101AN/BN: cheaper learning/prototype part, 100 W, 1.8-250 MHz, 50 V class part.
- NXP MRF300AN/BN: serious one-axis candidate, 300 W, 1.8-250 MHz, 50 V class part, mirror pinouts for push-pull.
- Bigger escalation candidates: MRFE6VP5600H, MRFX1K80H, BLF188XR/BLF578XR class parts. Do not start here unless needed.

Recommended prototype sequence:

1. One-axis dummy-load model in Jupyter.
2. One-axis dummy-load RF breadboard/PCB with capacitance matching measured crystal fixture.
3. One-axis optical test.
4. Duplicate for x/y quadrature.

