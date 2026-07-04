# RF Design Problem

## Top-level electrical model

Each EO axis is initially modeled as a differential capacitive load:

```math
Z_C(f)=\frac{1}{j2\pi f C_e}
```

Current requirement:

```math
I_{peak}=2\pi f C_e V_{diff,peak}
```

For a chirp, use instantaneous frequency or compute current from time-domain derivative:

```math
i(t)=C_e\frac{dV_{diff}(t)}{dt}
```

## Current numbers to keep in mind

With `C_e = 13 pF`:

- `Vdiff = 480 Vpeak`, `f = 10 MHz` -> `Ipeak ≈ 0.39 A`
- `Vdiff = 480 Vpeak`, `f = 50 MHz` -> `Ipeak ≈ 1.96 A`
- `Vdiff = 480 Vpeak`, `f = 100 MHz` -> `Ipeak ≈ 3.92 A`

This is the electrode current. If a transformer is used with voltage step-up `N`, the primary current is roughly `N * Is` for an ideal transformer.

## Why a big output step-up transformer is dangerous

For voltage ratio `N = Vs/Vp`:

```math
Z_p = \frac{Z_s}{N^2}, \quad C_{seen}=N^2 C_s
```

So a 1:3 voltage step-up makes a 13 pF load look like about 117 pF at the primary. This is why a big output transformer can make the amplifier current/reflection problem worse.

## Preferred PA design stance

Use an LDMOS/GaN RF PA where the EO capacitance is part of the final output network. Do not design a generic 50-ohm PA and then attach the EO capacitor as an afterthought.

For one axis:

```text
chirp source
  -> driver / phase splitter
  -> push-pull LDMOS PA
  -> output match/equalizer designed with EO capacitance included
  -> EO electrode pair
  -> capacitive RF pickup
```

## What the notebook should optimize

For each candidate output network, compute vs frequency:

- electrode voltage magnitude `|V_EO(f)|`
- electrode-voltage phase `angle(V_EO(f))`
- current through EO capacitance
- current through damping elements
- real power dissipated in damping/termination elements
- equivalent load seen by transistor/output stage
- smoothness / monotonicity of transfer function for chirp predistortion

The goal is not necessarily a perfect 50-ohm match. The goal is a high electrode voltage with manageable transistor current/voltage, manageable heating, and repeatable transfer function.

## Practical path

1. Build a dummy load approximating the measured EO capacitance and fixture parasitics.
2. Model it in Jupyter.
3. Test output networks in simulation:
   - naked capacitive load
   - series damping
   - shunt damping
   - low-Q peaking/equalization
   - modest transformer/balun only if needed
4. Measure with VNA / directional coupler / RF pickup.
5. Fit measured S-parameters to the model.
6. Apply predistortion to the AWG waveform.

