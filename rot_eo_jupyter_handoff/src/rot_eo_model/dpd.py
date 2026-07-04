"""Predistortion and transfer-function utilities.

This is an initial linear scaffold. Later versions can include nonlinear/ML DPD.
"""

from __future__ import annotations

import numpy as np


def safe_inverse_transfer(H: np.ndarray, floor: float = 1e-3) -> np.ndarray:
    """Return a regularized inverse transfer function 1/H.

    Parameters
    ----------
    H:
        Complex transfer function.
    floor:
        Minimum allowed magnitude to avoid exploding the inverse.
    """
    H = np.asarray(H, dtype=complex)
    mag = np.abs(H)
    H_reg = np.where(mag < floor, floor * np.exp(1j * np.angle(H)), H)
    return 1.0 / H_reg


def predistort_spectrum(target_spectrum: np.ndarray, H: np.ndarray, max_gain: float = 20.0) -> np.ndarray:
    """Predistort a target spectrum by inverse H with gain clipping."""
    inv = safe_inverse_transfer(H)
    gain = np.abs(inv)
    inv_clipped = np.where(gain > max_gain, inv * (max_gain / gain), inv)
    return target_spectrum * inv_clipped


def quadrature_error_metrics(Hx: np.ndarray, Hy: np.ndarray) -> dict[str, np.ndarray]:
    """Return amplitude ratio and phase error between x and y channels.

    Phase error is angle(Hy/Hx) - pi/2, wrapped to [-pi, pi].
    """
    Hx = np.asarray(Hx, dtype=complex)
    Hy = np.asarray(Hy, dtype=complex)
    amp_ratio = np.abs(Hy) / np.maximum(np.abs(Hx), 1e-30)
    phase_error = np.angle(Hy / Hx) - np.pi / 2.0
    phase_error = (phase_error + np.pi) % (2.0 * np.pi) - np.pi
    return {"amp_ratio_y_over_x": amp_ratio, "phase_error_rad": phase_error}
