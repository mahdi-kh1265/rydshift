# rot_eo_ads_spice_export_v1_500um

## What is this package?
This is a curated, versioned export of the RF-ready linear load models for the rotating electro-optic (EO) frequency shifter. It contains machine-readable Touchstone (.sNp) and SPICE (.sp/.sub) files.

## What geometry does it represent?
It represents the locked 500 µm hybrid setback geometry (z-cut LiNbO3, 3x3x30 mm).

## Which file should I use first in ADS?
Use `ADS/eo_load_4port_single_ended_R50.s4p` using the SnP component. See `ADS_FIRST_USE_CHECKLIST.md`.

## Which file should I use first in LTspice/ngspice?
Use `SPICE_LTspice_ngspice/eo_load_4node.sp`. See `SPICE_FIRST_USE_CHECKLIST.md`.

## Port Definitions
See `geometry_source/port_map.md` for definitions.

## Differences between models
- **Lossless**: Pure authoritative capacitance matrix.
- **Lossy (RC-Ladder)**: Passive broadband approximation of a constant loss tangent.
- **Tolerance**: Scaled bounds (±20%) for extreme operating corner checking.
- **Fixture Wrapper**: A template for appending physical parasitics (wirebonds/traces).

## Validation Status
- **Math & Formatting**: Validated perfectly against Touchstone/SPICE specifications.
- **ADS/LTspice Smoke Testing**: **PENDING** actual simulator import by the RF team.
