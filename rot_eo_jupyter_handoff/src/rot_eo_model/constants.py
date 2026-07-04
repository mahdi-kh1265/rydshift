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
    capacitance_F: float = 13e-12
    f_min_Hz: float = 6e6
    f_max_Hz: float = 102e6
    vpi_diff_peak_V: float = 480.0

DEFAULTS = EOConstants()
