"""Electro-optic tensor and retardance utilities for z-cut LiNbO3.

This module implements the simplified z-propagation transverse impermeability
model discussed in the handoff using rigorous numerical diagonalization.
"""

from __future__ import annotations

import numpy as np
from rot_eo_model.tensors import PockelsTensor


def eta_o(n_o: float) -> float:
    """Return ordinary impermeability eta_o = 1 / n_o^2."""
    return 1.0 / (n_o * n_o)


def transverse_eta_matrix(
    Ex_V_per_m: float,
    Ey_V_per_m: float,
    Ez_V_per_m: float = 0.0,
    *,
    n_o: float = 2.286,
    pockels: PockelsTensor,
) -> np.ndarray:
    """Return the 2x2 transverse impermeability matrix for z-propagation.

    Parameters
    ----------
    Ex_V_per_m, Ey_V_per_m, Ez_V_per_m:
        Electric-field components in crystal axes.
    n_o:
        Ordinary refractive index.
    pockels:
        Pockels tensor containing r_ijk.

    Returns
    -------
    eta_perp : ndarray, shape (2, 2)
        Symmetric transverse impermeability matrix.
    """
    e0 = eta_o(n_o)
    
    # Contracted indices: 0:xx, 1:yy, 2:zz, 3:yz, 4:xz, 5:xy
    r = pockels.contracted
    
    # d(eta_i) = r_ij E_j
    # We only care about eta_1 (xx), eta_2 (yy), and eta_6 (xy) for transverse z-propagation
    d_eta_xx = r[0, 0] * Ex_V_per_m + r[0, 1] * Ey_V_per_m + r[0, 2] * Ez_V_per_m
    d_eta_yy = r[1, 0] * Ex_V_per_m + r[1, 1] * Ey_V_per_m + r[1, 2] * Ez_V_per_m
    d_eta_xy = r[5, 0] * Ex_V_per_m + r[5, 1] * Ey_V_per_m + r[5, 2] * Ez_V_per_m
    
    return np.array(
        [
            [e0 + d_eta_xx, d_eta_xy],
            [d_eta_xy, e0 + d_eta_yy],
        ],
        dtype=float,
    )


def eigen_indices(eta_perp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return refractive indices and eigenvectors of a transverse eta matrix.

    Eigenvalues lambda_i are equal to 1/n_i^2.
    """
    vals, vecs = np.linalg.eigh(np.asarray(eta_perp, dtype=float))
    if np.any(vals <= 0):
        raise ValueError("Impermeability eigenvalues must be positive.")
    n = 1.0 / np.sqrt(vals)
    return n, vecs


def full_retardance_eigenvalue(
    E_perp_V_per_m: np.ndarray,
    *,
    length_m: float = 30e-3,
    lambda0_m: float = 780e-9,
    n_o: float = 2.286,
    pockels: PockelsTensor,
) -> float:
    """Exact retardance calculation from numerical eigenvalues.
    
    E_perp_V_per_m is [Ex, Ey, Ez].
    """
    eta = transverse_eta_matrix(
        E_perp_V_per_m[0], E_perp_V_per_m[1], E_perp_V_per_m[2],
        n_o=n_o, pockels=pockels
    )
    n, _ = eigen_indices(eta)
    dn = abs(n[1] - n[0])
    return 2.0 * np.pi / lambda0_m * length_m * dn


def delta_n_small_signal(
    E_perp_V_per_m: float,
    *,
    n_o: float = 2.286,
    r22_m_per_V: float = 6.8e-12,
) -> float:
    """Approximate birefringent splitting Δn for z-cut rotating EO geometry."""
    return n_o**3 * abs(r22_m_per_V) * abs(E_perp_V_per_m)


def retardance_rad(
    E_perp_V_per_m: float,
    *,
    lambda0_m: float = 780e-9,
    length_m: float = 30e-3,
    n_o: float = 2.286,
    r22_m_per_V: float = 6.8e-12,
) -> float:
    """Return small-signal retardance Γ in radians."""
    dn = delta_n_small_signal(E_perp_V_per_m, n_o=n_o, r22_m_per_V=r22_m_per_V)
    return 2.0 * np.pi / lambda0_m * length_m * dn


def half_wave_field_V_per_m(
    *,
    lambda0_m: float = 780e-9,
    length_m: float = 30e-3,
    n_o: float = 2.286,
    r22_m_per_V: float = 6.8e-12,
) -> float:
    """Return transverse field E_pi required for Γ = pi."""
    return lambda0_m / (2.0 * length_m * n_o**3 * abs(r22_m_per_V))


def vpi_from_gap_V(
    gap_m: float,
    *,
    lambda0_m: float = 780e-9,
    length_m: float = 30e-3,
    n_o: float = 2.286,
    r22_m_per_V: float = 6.8e-12,
    field_efficiency: float = 1.0,
) -> float:
    """Return differential Vpi using E≈field_efficiency*V/gap."""
    Epi = half_wave_field_V_per_m(
        lambda0_m=lambda0_m, length_m=length_m, n_o=n_o, r22_m_per_V=r22_m_per_V
    )
    return Epi * gap_m / field_efficiency


def axis_angle_rad(vecs: np.ndarray) -> float:
    """Return principal-axis angle theta from eigenvectors.
    
    Extracts the angle of the principal axis (first eigenvector) 
    with respect to the x-axis, using exact numerical vectors.
    """
    v1 = vecs[:, 0]
    return np.arctan2(v1[1], v1[0])
