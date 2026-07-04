"""Project constants for rotating EO LiNbO3 model.

All values are SI unless explicitly stated.
These defaults are working assumptions, not final measured truth.
"""

from dataclasses import dataclass
from rot_eo_model.materials import linbo3_pockels_tensor, linbo3_dielectric_tensor, sellmeier_no
from rot_eo_model.tensors import PockelsTensor, DielectricTensor

@dataclass(frozen=True)
class EOConstants:
    lambda0_m: float = 780e-9
    
    @property
    def n_o(self) -> float:
        return sellmeier_no(self.lambda0_m)
    
    r22_m_per_V: float = 6.8e-12
    r13_m_per_V: float = 9.6e-12
    
    @property
    def pockels(self) -> PockelsTensor:
        return linbo3_pockels_tensor(r22_m_per_V=self.r22_m_per_V, r13_m_per_V=self.r13_m_per_V)
    
    crystal_length_m: float = 30e-3
    electrode_gap_m: float = 3e-3
    
    # Locked hybrid geometry: 500 um setback, full-face +x/-x, 2 mm strip +y/-y
    strip_setback_m: float = 500e-6
    strip_width_m: float = 2e-3
    
    # FEM-derived differential capacitances (asymmetric channels)
    capacitance_xx_F: float = 16.880e-12   # full-face channel
    capacitance_yy_F: float = 16.179e-12   # strip channel
    capacitance_xy_F: float = 0.0          # negligible cross-coupling
    
    # Backward compatibility: use the larger (bottleneck) channel
    capacitance_F: float = 16.880e-12
    
    # FEM-derived half-wave voltages (asymmetric)
    vpi_x_V: float = 571.9
    vpi_y_V: float = 628.5
    vpi_diff_peak_V: float = 628.5  # use the bottleneck (y-channel)
    
    f_min_Hz: float = 6e6
    f_max_Hz: float = 102e6

DEFAULTS = EOConstants()
