"""Generate ADS and SPICE exports for the EO load from the JSON specification."""

import json
import numpy as np
import skrf as rf
from scipy.optimize import nnls
from pathlib import Path
import matplotlib.pyplot as plt

def fit_rc_ladder(freqs, tand, n_poles=6):
    """
    Fit a normalized passive RC-ladder to approximate a constant tan(delta)
    capacitor over the given frequencies. Target C_norm = 1.0 Farads.
    """
    omega = 2 * np.pi * freqs
    w_min = omega[0] / 3.0
    w_max = omega[-1] * 3.0
    
    omega_poles = np.logspace(np.log10(w_min), np.log10(w_max), n_poles)
    tau = 1.0 / omega_poles
    
    M = len(omega)
    A = np.zeros((2*M, n_poles + 1))
    b = np.zeros(2*M)
    
    for i, w in enumerate(omega):
        A[i, 0] = 1.0
        for k in range(n_poles):
            A[i, k+1] = 1.0 / (1.0 + (w * tau[k])**2)
        b[i] = 1.0
        
        A[M+i, 0] = 0.0
        for k in range(n_poles):
            A[M+i, k+1] = (w * tau[k]) / (1.0 + (w * tau[k])**2)
        weight = 1.0 / tand 
        A[M+i, :] *= weight
        b[M+i] = tand * weight

    x, _ = nnls(A, b)
    C_inf = x[0]
    C_k = x[1:]
    R_k = tau / C_k
    return C_inf, C_k, R_k

def main():
    root_dir = Path(__file__).parent.parent
    ads_dir = Path(__file__).parent
    
    # 1. Load spec
    spec_path = root_dir / "rf_ready_geometry_spec.json"
    with open(spec_path, "r") as f:
        spec = json.load(f)
        
    c_node = np.array(spec["C_node_pF"]) * 1e-12
    c_diff = np.array(spec["C_diff_pF"]) * 1e-12
    
    # Validation checks
    assert np.allclose(c_node, c_node.T), "C_node must be symmetric"
    assert np.allclose(c_diff, c_diff.T), "C_diff must be symmetric"
    
    # =========================================================================
    # Generate SPICE Subcircuits
    # =========================================================================
    def write_spice_subckt(filename, c_mat, tand=0.0):
        with open(ads_dir / filename, "w") as f:
            f.write(f"* EO Load Subcircuit: {filename}\n")
            f.write(f"* Loss tangent: {tand}\n")
            if tand > 0:
                f.write("* WARNING: This is a SINGLE-FREQUENCY approximation evaluated exactly at 100 MHz.\n")
                f.write("* It does NOT represent a constant tan(delta) over frequency.\n")
            f.write(f"* Original Maxwell C_node matrix (F):\n")
            for row in c_mat:
                f.write("* " + "  ".join(f"{v:.4e}" for v in row) + "\n")
            f.write("\n")
            f.write(".SUBCKT EO_LOAD_4NODE p_x m_x p_y m_y\n")
            
            nodes = ["p_x", "m_x", "p_y", "m_y"]
            n = 4
            idx = 1
            
            for i in range(n):
                for j in range(i+1, n):
                    c_val = -c_mat[i, j]
                    if c_val > 1e-18:
                        f.write(f"C{idx} {nodes[i]} {nodes[j]} {c_val:.6e}\n")
                        if tand > 0:
                            r_val = 1.0 / (2 * np.pi * 100e6 * c_val * tand)
                            f.write(f"R{idx} {nodes[i]} {nodes[j]} {r_val:.6e}\n")
                        idx += 1
                        
            for i in range(n):
                c_gnd = np.sum(c_mat[i, :])
                if c_gnd > 1e-18:
                    f.write(f"C{idx} {nodes[i]} 0 {c_gnd:.6e}\n")
                    if tand > 0:
                        r_val = 1.0 / (2 * np.pi * 100e6 * c_gnd * tand)
                        f.write(f"R{idx} {nodes[i]} 0 {r_val:.6e}\n")
                    idx += 1
            f.write(".ENDS\n")

    def write_rc_ladder_subckt(filename, c_mat, tand):
        freqs = np.linspace(0.5e6, 200e6, 400)
        C_inf, C_k, R_k = fit_rc_ladder(freqs, tand, n_poles=6)
        
        with open(ads_dir / filename, "w") as f:
            f.write(f"* Broadband RC-Ladder EO Load Subcircuit: {filename}\n")
            f.write(f"* Loss tangent: {tand}\n")
            f.write("* This model approximates constant tan(delta) over 0.5-200 MHz.\n")
            f.write("* It is a finite-band approximation, not an exact all-frequency constant-tand law.\n\n")
            f.write(".SUBCKT EO_LOAD_4NODE p_x m_x p_y m_y\n")
            
            nodes = ["p_x", "m_x", "p_y", "m_y"]
            n = 4
            idx = 1
            
            def write_ladder(node1, node2, c_val):
                nonlocal idx
                f.write(f"* Edge: {node1} to {node2}, base C = {c_val:.4e}\n")
                f.write(f"Cinf_{idx} {node1} {node2} {c_val * C_inf:.6e}\n")
                for k in range(len(C_k)):
                    if C_k[k] > 1e-15:
                        f.write(f"Ck_{idx}_{k} {node1} n_int_{idx}_{k} {c_val * C_k[k]:.6e}\n")
                        # Resistance is inversely proportional to scaling
                        f.write(f"Rk_{idx}_{k} n_int_{idx}_{k} {node2} {R_k[k] / c_val:.6e}\n")
                idx += 1

            for i in range(n):
                for j in range(i+1, n):
                    c_val = -c_mat[i, j]
                    if c_val > 1e-18:
                        write_ladder(nodes[i], nodes[j], c_val)
                        
            for i in range(n):
                c_gnd = np.sum(c_mat[i, :])
                if c_gnd > 1e-18:
                    write_ladder(nodes[i], "0", c_gnd)
                    
            f.write(".ENDS\n")

    # Generate purely capacitive SPICE model
    write_spice_subckt("eo_load_4node.sp", c_node, tand=0.0)
    
    # Generate 100MHz-only lossy variants and broadband variants
    for td in [1e-4, 5e-4, 1e-3]:
        write_spice_subckt(f"eo_load_4node_lossy_tand_{td:g}_at_100MHz.sp", c_node, tand=td)
        write_rc_ladder_subckt(f"eo_load_4node_lossy_tand_{td:g}_ladder_0p5MHz_200MHz.sub", c_node, tand=td)

    # =========================================================================
    # Generate Touchstone Files
    # =========================================================================
    freqs = np.linspace(0.5e6, 200e6, 400)
    freq = rf.Frequency.from_f(freqs, unit='hz')
    omega = 2 * np.pi * freqs[:, np.newaxis, np.newaxis]
    
    def make_network(name, c_mat, z0, tand=0.0):
        Y = (omega * tand + 1j * omega) * c_mat[np.newaxis, :, :]
        ntwk = rf.Network(frequency=freq, y=Y, name=name)
        ntwk.renormalize(z0)
        
        # Validation
        assert np.allclose(ntwk.s, np.transpose(ntwk.s, axes=(0, 2, 1))), "S-matrix must be reciprocal"
        if tand >= 0:
            s_svd = np.linalg.svd(ntwk.s, compute_uv=False)
            assert np.max(s_svd) <= 1.0001, "Network must be passive"
            
        return ntwk

    # 4-port single-ended (50 ohm)
    ntwk_4p = make_network("eo_load_4port_single_ended_R50", c_node, z0=50.0, tand=0.0)
    ntwk_4p.write_touchstone(filename="eo_load_4port_single_ended_R50.s4p", dir=str(ads_dir))
    
    # Tolerances
    for scale in [0.8, 1.2]:
        ntwk_4p_tol = make_network(f"eo_load_4port_single_ended_R50_tol_{scale}", c_node * scale, z0=50.0, tand=0.0)
        ntwk_4p_tol.write_touchstone(filename=f"eo_load_4port_single_ended_R50_tol_{scale}.s4p", dir=str(ads_dir))
    
    for td in [1e-4, 5e-4, 1e-3]:
        ntwk_4p_lossy = make_network(f"eo_load_4port_single_ended_R50_tand_{td:g}", c_node, z0=50.0, tand=td)
        ntwk_4p_lossy.write_touchstone(filename=f"eo_load_4port_single_ended_R50_tand_{td:g}.s4p", dir=str(ads_dir))
        
    # 2-port differential
    ntwk_2p_50 = make_network("eo_load_2port_diff_R50", c_diff, z0=50.0, tand=0.0)
    ntwk_2p_50.write_touchstone(filename="eo_load_2port_diff_R50.s2p", dir=str(ads_dir))
    
    ntwk_2p_100 = make_network("eo_load_2port_diff_R100", c_diff, z0=100.0, tand=0.0)
    ntwk_2p_100.write_touchstone(filename="eo_load_2port_diff_R100.s2p", dir=str(ads_dir))
    
    for td in [1e-4, 5e-4, 1e-3]:
        ntwk_2p_lossy = make_network(f"eo_load_2port_diff_R50_tand_{td:g}", c_diff, z0=50.0, tand=td)
        ntwk_2p_lossy.write_touchstone(filename=f"eo_load_2port_diff_R50_tand_{td:g}.s2p", dir=str(ads_dir))
        
    # 1-port differential for x and y independently
    cx = c_diff[0:1, 0:1]
    cy = c_diff[1:2, 1:2]
    make_network("eo_load_xdiff_R50", cx, z0=50.0).write_touchstone(filename="eo_load_xdiff_R50.s1p", dir=str(ads_dir))
    make_network("eo_load_ydiff_R50", cy, z0=50.0).write_touchstone(filename="eo_load_ydiff_R50.s1p", dir=str(ads_dir))

    # =========================================================================
    # Validation Plots
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot effective tan(delta) of RC ladders
    ax = axes[0, 0]
    for td in [1e-4, 5e-4, 1e-3]:
        C_inf, C_k, R_k = fit_rc_ladder(freqs, td, n_poles=6)
        Y_norm = 1j * omega.flatten() * C_inf
        for k in range(len(C_k)):
            if C_k[k] > 1e-15:
                Y_norm += (1j * omega.flatten() * C_k[k]) / (1 + 1j * omega.flatten() * R_k[k] * C_k[k])
        tand_eff = np.real(Y_norm) / np.imag(Y_norm)
        ax.plot(freqs/1e6, tand_eff, label=f'Ladder (target {td:g})')
        ax.axhline(td, color='k', linestyle='--', alpha=0.3)
    ax.set_title("RC-Ladder Effective tan(δ)")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("tan(δ)")
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True)
    
    # Tolerance impedance magnitude
    ax = axes[0, 1]
    n_08 = make_network("n_08", c_diff * 0.8, 50.0)
    n_12 = make_network("n_12", c_diff * 1.2, 50.0)
    
    ax.plot(freq.f/1e6, np.abs(ntwk_2p_50.z[:,0,0]), label='Nominal x diff')
    ax.plot(freq.f/1e6, np.abs(n_08.z[:,0,0]), '--', label='0.8x tol')
    ax.plot(freq.f/1e6, np.abs(n_12.z[:,0,0]), '--', label='1.2x tol')
    ax.set_title("Tolerance Bounds on |Z|")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Ohms")
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True)
    
    # Current at Vpi
    ax = axes[1, 0]
    Vpi_x = spec["Vpi_x_V"]
    ax.plot(freq.f/1e6, Vpi_x / np.abs(ntwk_2p_50.z[:,0,0]), label='Nominal Ix')
    ax.plot(freq.f/1e6, Vpi_x / np.abs(n_08.z[:,0,0]), '--', label='0.8x Ix')
    ax.plot(freq.f/1e6, Vpi_x / np.abs(n_12.z[:,0,0]), '--', label='1.2x Ix')
    ax.set_title("Peak Current at $V_\\pi$ (x diff)")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Amps")
    ax.legend()
    ax.grid(True)
    
    # Ladder Z Phase
    ax = axes[1, 1]
    ax.plot(freq.f/1e6, np.angle(ntwk_2p_50.z[:,0,0], deg=True), label='Lossless')
    # Phase for lossy diff
    td = 1e-3
    n_lossy = rf.Network(str(ads_dir / f"eo_load_2port_diff_R50_tand_{td:g}.s2p"))
    ax.plot(freq.f/1e6, np.angle(n_lossy.z[:,0,0], deg=True), '--', label=f'Lossy tand={td}')
    ax.set_title("Input Impedance Phase (x diff)")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Degrees")
    ax.set_ylim(-91, -89)
    ax.legend()
    ax.grid(True)
    
    fig.tight_layout()
    fig.savefig(ads_dir / "validation_plots.png", dpi=150)
    
    # =========================================================================
    # Write README
    # =========================================================================
    readme = r"""# ADS / SPICE Import Package for Rotating EO Modulator

This directory contains the RF-ready linear load models for the locked 500 µm hybrid geometry.

> **WARNING**: These are strictly small-signal linear load models. They do **NOT** model the nonlinear LDMOS PA or the matching network.

## Touchstone Files (.sNp)
Use these in ADS using the `SnP` component. The files are written in **RI** (Real/Imaginary) format.
Includes `_tol_0.8` and `_tol_1.2` tolerance variants for the 4-port core.

## Minimal ADS SnP Smoke Test
To verify format compliance in ADS before proceeding to layout/matching:
1. Place an `SnP` component in an empty ADS schematic.
2. Point it to `eo_load_4port_single_ended_R50.s4p`.
3. Terminate all four ports with 50 Ohm (`Term` components).
4. Run an S-parameter sweep from 0.5 MHz to 200 MHz.
5. Check expected input capacitance (Y-parameters) at low frequency:
   - $C_{diff,x} \approx 16.88$ pF
   - $C_{diff,y} \approx 16.18$ pF
6. Check input impedance magnitude $|Z|$ at 100 MHz:
   - $|Z_{diff,x}| \approx 94\ \Omega$
   - $|Z_{diff,y}| \approx 98\ \Omega$

**IMPORTANT**: The 4-port `.s4p` is the definitive authoritative physical model. All other `.s2p` models are derivative reductions.

## SPICE Subcircuits
`eo_load_4node.sp` is the authoritative pure-capacitance Maxwell model.

**Important Note on Lossy Models**:
- The `_at_100MHz.sp` files contain fixed resistors valid EXACTLY at 100 MHz.
- The `_ladder_0p5MHz_200MHz.sub` files are broadband passive RC-ladder approximations that maintain a nearly constant $\tan\delta$ over the specified bandwidth. These are robust for transient and AC analysis.

## Fixture Parasitics Wrapper
The `fixture_parasitics_wrapper.net` (ADS syntax) demonstrates how to properly wrap the ideal 4-port core with estimated series inductances and resistances. *These parasitics must be verified against VNA measurements.*

## Measurement Fitting
`fit_vna_measurements.py` provides a scaffold to extract $C_{node}$, $\tan\delta$, and series $R/L$ directly from VNA `.s4p` files.
"""
    with open(ads_dir / "README_ADS_IMPORT.md", "w", encoding="utf-8") as f:
        f.write(readme)
        
    print("Export complete.")

if __name__ == "__main__":
    main()
