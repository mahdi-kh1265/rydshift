"""
Measurement Fitting Scaffold

Reads measured VNA .s4p files, converts to Y-parameters, and extracts:
- Effective C_node Maxwell capacitance matrix
- Effective loss tangent
- Series fixture parasitics (R_s, L_s)

This script is meant to be run by the measurement/RF team when actual VNA data is available.
"""
import numpy as np
import skrf as rf
from pathlib import Path
import json

def fit_low_frequency_capacitance(ntwk, f_max_Hz=50e6):
    """
    Extract the effective 4x4 C_node matrix from the low-frequency imaginary Y-parameters.
    For Y = j*w*C, C = Im(Y) / w.
    """
    f_mask = ntwk.f <= f_max_Hz
    y_low = ntwk.y[f_mask]
    omega = 2 * np.pi * ntwk.f[f_mask][:, np.newaxis, np.newaxis]
    
    # C(w) = Im(Y(w)) / w
    C_w = np.imag(y_low) / omega
    
    # Average over the low-frequency band
    C_fit = np.mean(C_w, axis=0)
    
    # Enforce symmetry
    C_fit = 0.5 * (C_fit + C_fit.T)
    return C_fit

def fit_loss_tangent(ntwk, f_min_Hz=10e6, f_max_Hz=100e6):
    """
    Extract effective dielectric loss tangent from Re(Y) / Im(Y).
    """
    f_mask = (ntwk.f >= f_min_Hz) & (ntwk.f <= f_max_Hz)
    y_mid = ntwk.y[f_mask]
    
    # Avoid dividing by zero or extracting from cross-coupling terms that are ~0
    tand_w = np.real(y_mid) / np.abs(np.imag(y_mid) + 1e-15)
    
    # Only fit the diagonal terms where the signal is strong
    tand_diag = np.diagonal(tand_w, axis1=1, axis2=2)
    tand_fit = np.mean(tand_diag, axis=0)  # average per port
    return np.mean(tand_fit)  # global average

def fit_series_parasitics(ntwk, C_ideal, f_min_Hz=150e6):
    """
    At high frequencies, the series parasitics dominate.
    Z_meas = Z_ideal + Z_series = Z_ideal + R_s + j*w*L_s
    Z_series = Z_meas - inv(Y_ideal)
    """
    f_mask = ntwk.f >= f_min_Hz
    z_high = ntwk.z[f_mask]
    omega = 2 * np.pi * ntwk.f[f_mask][:, np.newaxis, np.newaxis]
    
    # Ideal admittance of pure capacitance
    Y_ideal = 1j * omega * C_ideal[np.newaxis, :, :]
    
    # Z_ideal = inv(Y_ideal)
    Z_ideal = np.linalg.inv(Y_ideal)
    
    Z_series = z_high - Z_ideal
    
    R_s = np.mean(np.real(Z_series), axis=0)
    L_s = np.mean(np.imag(Z_series) / omega, axis=0)
    
    # Only return diagonal (per-port) parasitics
    return np.diag(R_s), np.diag(L_s)

def main():
    print("--- VNA Measurement Fitting Scaffold ---")
    
    # Load the nominal .s4p as a dummy "measured" file for demonstration
    base_dir = Path(__file__).parent
    dummy_file = base_dir / "eo_load_4port_single_ended_R50.s4p"
    
    if not dummy_file.exists():
        print(f"Error: {dummy_file} not found.")
        return
        
    ntwk = rf.Network(str(dummy_file))
    print(f"Loaded {ntwk.name} ({ntwk.f[0]/1e6:.1f} - {ntwk.f[-1]/1e6:.1f} MHz)")
    
    # 1. Fit C_node
    C_fit = fit_low_frequency_capacitance(ntwk)
    print("\nFitted C_node matrix (pF):")
    print(np.round(C_fit * 1e12, 3))
    
    # Load FEM baseline
    spec_path = base_dir.parent / "rf_ready_geometry_spec.json"
    if spec_path.exists():
        with open(spec_path, "r") as f:
            spec = json.load(f)
        C_fem = np.array(spec["C_node_pF"]) * 1e-12
        err = np.max(np.abs(C_fit - C_fem))
        print(f"Max error vs FEM baseline: {err * 1e12:.5f} pF")
        
    # 2. Fit loss tangent
    tand_fit = fit_loss_tangent(ntwk)
    print(f"\nFitted loss tangent: {tand_fit:.2e}")
    
    # 3. Fit series parasitics
    # (Since this is a dummy lossless file with no parasitics, these should be ~0)
    Rs, Ls = fit_series_parasitics(ntwk, C_fit)
    print(f"\nFitted Series R (Ohms): {np.round(Rs, 3)}")
    print(f"Fitted Series L (nH):   {np.round(Ls * 1e9, 3)}")
    
    print("\nScaffold complete. Ready for real VNA data.")

if __name__ == "__main__":
    main()
