# SPICE / LTspice / ngspice First Use Checklist

- **Pure capacitor model** (`eo_load_4node.sp`) is the portable baseline.
- **Fixed resistor lossy models** (`...at_100MHz.sp`) are exact **only** at 100 MHz.
- **RC-ladder lossy models** (`...ladder_0p5MHz_200MHz.sub`) are finite-band approximations over 0.5–200 MHz, great for transient sweeps.
- To instantiate: Use `.INCLUDE "eo_load_4node.sp"` and instance it as `X1 p_x m_x p_y m_y EO_LOAD_4NODE`.
- **LTspice Status**: Pending actual simulator import.
