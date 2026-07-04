"""Lumped RF utilities for EO capacitive-load screening.

These are simple, transparent calculations intended to set specs and kill
bad architectures early. They are not a replacement for nonlinear RF PA design.
"""

from __future__ import annotations

import numpy as np


def z_capacitor(freq_Hz: np.ndarray | float, C_F: float) -> np.ndarray:
    """Return capacitor impedance 1/(jωC)."""
    f = np.asarray(freq_Hz, dtype=float)
    return 1.0 / (1j * 2.0 * np.pi * f * C_F)


def cap_current_peak_A(freq_Hz: np.ndarray | float, C_F: float, Vpeak_V: float) -> np.ndarray:
    """Peak sinusoidal current through a capacitor."""
    f = np.asarray(freq_Hz, dtype=float)
    return 2.0 * np.pi * f * C_F * Vpeak_V


def cap_apparent_power_VA(freq_Hz: np.ndarray | float, C_F: float, Vpeak_V: float) -> np.ndarray:
    """Reactive apparent power scale Vrms*Irms for sinusoidal capacitor drive."""
    Ipk = cap_current_peak_A(freq_Hz, C_F, Vpeak_V)
    return (Vpeak_V / np.sqrt(2.0)) * (Ipk / np.sqrt(2.0))


def dielectric_loss_W(freq_Hz: np.ndarray | float, C_F: float, Vpeak_V: float, tan_delta: float) -> np.ndarray:
    """Approximate real dielectric loss in a capacitor with loss tangent tanδ.

    P = 0.5 * ω C V_peak^2 tanδ
    """
    f = np.asarray(freq_Hz, dtype=float)
    return 0.5 * 2.0 * np.pi * f * C_F * Vpeak_V**2 * tan_delta


def reflected_capacitance_F(C_secondary_F: float, turns_ratio_voltage: float) -> float:
    """Effective primary capacitance for an ideal step-up transformer.

    If Vs = N*Vp, then impedance reflects as Zp=Zs/N^2 and capacitance
    reflects as Cp=N^2*Cs.
    """
    N = turns_ratio_voltage
    return N * N * C_secondary_F


def reflected_impedance(Z_secondary: np.ndarray | complex, turns_ratio_voltage: float) -> np.ndarray:
    """Reflect secondary impedance to primary through ideal transformer."""
    return np.asarray(Z_secondary) / (turns_ratio_voltage**2)


def sinusoidal_waveform_current(t_s: np.ndarray, Vdiff_V: np.ndarray, C_F: float) -> np.ndarray:
    """Time-domain capacitor current i=C*dV/dt using numerical derivative."""
    t = np.asarray(t_s, dtype=float)
    v = np.asarray(Vdiff_V, dtype=float)
    dvdt = np.gradient(v, t)
    return C_F * dvdt


def shunt_parallel(Z1: complex | np.ndarray, Z2: complex | np.ndarray) -> np.ndarray:
    """Parallel combination of two impedances."""
    return 1.0 / (1.0 / np.asarray(Z1) + 1.0 / np.asarray(Z2))


def series(*Zs: complex | np.ndarray) -> np.ndarray:
    """Series sum of impedances."""
    total = 0
    for z in Zs:
        total = total + np.asarray(z)
    return total


def voltage_divider_load_voltage(Vsource_V: float, Zsource: complex | np.ndarray, Zload: complex | np.ndarray) -> np.ndarray:
    """Load voltage for Thevenin source with source impedance."""
    return Vsource_V * np.asarray(Zload) / (np.asarray(Zsource) + np.asarray(Zload))
