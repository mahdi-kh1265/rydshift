import numpy as np
from rot_eo_model.tensors import PockelsTensor, DielectricTensor

def sellmeier_no(lambda_m: float) -> float:
    """Ordinary refractive index of congruent LiNbO3 using Sellmeier equations.
    
    Formula based on Zelmon et al. (1997) or similar standard fits.
    At 780 nm, this should yield ~2.286.
    
    Wavelength lambda_m is in meters.
    """
    wl_um = lambda_m * 1e6
    # Jundt 1997 or Zelmon 1997 Sellmeier for congruent LN at 20 C
    # Using generalized form: n^2 = A + (B * l^2) / (l^2 - C) + (D * l^2) / (l^2 - E) - F * l^2
    # For simplicity, using a standard approximation for n_o that works well in NIR:
    # A = 4.9048, B = 0.11768, C = 0.04750, D = 0.027169
    # n_o^2 = 4.9048 + 0.11768 / (wl_um^2 - 0.04750) - 0.027169 * wl_um^2
    # Let's use the precise Jundt (1997) for congruent LN (at room temp ~ 24 C):
    # n_o^2 = 4.9048 + 0.11768 / (lambda_um^2 - 0.04750) - 0.027169 * lambda_um^2
    
    n_o_sq = 4.9048 + 0.11768 / (wl_um**2 - 0.04750) - 0.027169 * wl_um**2
    return np.sqrt(n_o_sq)

def sellmeier_ne(lambda_m: float) -> float:
    """Extraordinary refractive index of congruent LiNbO3 using Sellmeier equations.
    
    Wavelength lambda_m is in meters.
    """
    wl_um = lambda_m * 1e6
    # Jundt (1997) for congruent LN, n_e at room temp:
    # n_e^2 = 4.5820 + 0.099169 / (lambda_um^2 - 0.04443) - 0.02195 * lambda_um^2
    n_e_sq = 4.5820 + 0.099169 / (wl_um**2 - 0.04443) - 0.02195 * wl_um**2
    return np.sqrt(n_e_sq)

def linbo3_pockels_tensor(r22_m_per_V: float = 6.8e-12, r13_m_per_V: float = 9.6e-12) -> PockelsTensor:
    """Returns the contracted Pockels tensor (6x3) for point group 3m (LiNbO3).
    
    Usually for LN (point group 3m):
    r11 = 0, r12 = 0, r13 = r13
    r21 = 0, r22 = -r22, r23 = r13
    r31 = 0, r32 = 0, r33 = r33
    r41 = 0, r42 = r51, r43 = 0
    r51 = r51, r52 = 0, r53 = 0
    r61 = -r22, r62 = 0, r63 = 0
    
    Using the standard convention where the transverse terms are dominated by r22.
    """
    r33 = 30.9e-12
    r51 = 32.6e-12
    
    contracted = np.zeros((6, 3))
    # 0: 11
    contracted[0, 1] = -r22_m_per_V
    contracted[0, 2] = r13_m_per_V
    # 1: 22
    contracted[1, 1] = r22_m_per_V
    contracted[1, 2] = r13_m_per_V
    # 2: 33
    contracted[2, 2] = r33
    # 3: 23
    contracted[3, 1] = r51
    # 4: 13
    contracted[4, 0] = r51
    # 5: 12
    contracted[5, 0] = -r22_m_per_V
    
    return PockelsTensor(contracted=contracted)

def linbo3_dielectric_tensor(freq_Hz: float) -> DielectricTensor:
    """RF dielectric tensor for LiNbO3.
    
    At RF frequencies (6-100 MHz), the clamped/unclamped permittivity is roughly:
    epsilon_11 = epsilon_22 ≈ 44
    epsilon_33 ≈ 28
    """
    eps11 = 44.0
    eps22 = 44.0
    eps33 = 28.0
    return DielectricTensor.from_diagonal(eps11, eps22, eps33)

def linbo3_tan_delta(freq_Hz: float) -> float:
    """Approximate loss tangent (tan delta) of LiNbO3 at RF frequencies.
    
    Conservatively bound for RF design.
    """
    return 1e-3
