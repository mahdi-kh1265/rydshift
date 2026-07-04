"""Generate ADS and SPICE exports for the EO load from the JSON specification."""

import json
import numpy as np
import skrf as rf
from pathlib import Path
import matplotlib.pyplot as plt

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
        # Convert Maxwell C_node to pairwise components
        # C_ij = -C_node[i,j] for i != j
        # C_i0 = sum_j C_node[i,j]
        
        with open(ads_dir / filename, "w") as f:
            f.write(f"* EO Load Subcircuit: {filename}\n")
            f.write(f"* Loss tangent: {tand}\n")
            f.write(f"* Original Maxwell C_node matrix (F):\n")
            for row in c_mat:
                f.write("* " + "  ".join(f"{v:.4e}" for v in row) + "\n")
            f.write("\n")
            f.write(".SUBCKT EO_LOAD_4NODE p_x m_x p_y m_y\n")
            
            nodes = ["p_x", "m_x", "p_y", "m_y"]
            n = 4
            idx = 1
            
            # Pairwise caps
            for i in range(n):
                for j in range(i+1, n):
                    c_val = -c_mat[i, j]
                    if c_val > 1e-18:
                        f.write(f"C{idx} {nodes[i]} {nodes[j]} {c_val:.6e}\n")
                        if tand > 0:
                            # Equivalent parallel resistor at 100 MHz for simple SPICE compatibility
                            r_val = 1.0 / (2 * np.pi * 100e6 * c_val * tand)
                            f.write(f"R{idx} {nodes[i]} {nodes[j]} {r_val:.6e} * evaluated at 100 MHz for tand={tand}\n")
                        idx += 1
                        
            # Capacitance to ground
            for i in range(n):
                c_gnd = np.sum(c_mat[i, :])
                if c_gnd > 1e-18:
                    f.write(f"C{idx} {nodes[i]} 0 {c_gnd:.6e}\n")
                    if tand > 0:
                        r_val = 1.0 / (2 * np.pi * 100e6 * c_gnd * tand)
                        f.write(f"R{idx} {nodes[i]} 0 {r_val:.6e} * evaluated at 100 MHz for tand={tand}\n")
                    idx += 1
                    
            f.write(".ENDS\n")

    write_spice_subckt("eo_load_4node.sp", c_node, tand=0.0)
    for td in [1e-4, 5e-4, 1e-3]:
        write_spice_subckt(f"eo_load_4node_lossy_tand_{td:g}.sp", c_node, tand=td)

    # =========================================================================
    # Generate Touchstone Files
    # =========================================================================
    freqs = np.linspace(0.5e6, 200e6, 400)
    freq = rf.Frequency.from_f(freqs, unit='hz')
    omega = 2 * np.pi * freqs[:, np.newaxis, np.newaxis]
    
    def make_network(name, c_mat, z0, tand=0.0):
        # Y = (omega * tand + j*omega) * C
        Y = (omega * tand + 1j * omega) * c_mat[np.newaxis, :, :]
        ntwk = rf.Network(frequency=freq, y=Y, name=name)
        ntwk.z0 = z0
        
        # Validation
        assert np.allclose(ntwk.s, np.transpose(ntwk.s, axes=(0, 2, 1))), "S-matrix must be reciprocal"
        if tand >= 0:
            # Passivity: max singular value of S <= 1
            s_svd = np.linalg.svd(ntwk.s, compute_uv=False)
            assert np.max(s_svd) <= 1.0001, "Network must be passive"
            
        return ntwk

    # 4-port single-ended (50 ohm)
    ntwk_4p = make_network("eo_load_4port_single_ended_R50", c_node, z0=50.0, tand=0.0)
    ntwk_4p.write_touchstone(dir=str(ads_dir))
    
    for td in [1e-4, 5e-4, 1e-3]:
        ntwk_4p_lossy = make_network(f"eo_load_4port_single_ended_R50_tand_{td:g}", c_node, z0=50.0, tand=td)
        ntwk_4p_lossy.write_touchstone(dir=str(ads_dir))
        
    # 2-port differential (50 ohm)
    ntwk_2p_50 = make_network("eo_load_2port_diff_R50", c_diff, z0=50.0, tand=0.0)
    ntwk_2p_50.write_touchstone(dir=str(ads_dir))
    
    # 2-port differential (100 ohm)
    ntwk_2p_100 = make_network("eo_load_2port_diff_R100", c_diff, z0=100.0, tand=0.0)
    ntwk_2p_100.write_touchstone(dir=str(ads_dir))
    
    for td in [1e-4, 5e-4, 1e-3]:
        ntwk_2p_lossy = make_network(f"eo_load_2port_diff_R50_tand_{td:g}", c_diff, z0=50.0, tand=td)
        ntwk_2p_lossy.write_touchstone(dir=str(ads_dir))
        
    # 1-port differential for x and y independently
    cx = c_diff[0:1, 0:1]
    cy = c_diff[1:2, 1:2]
    make_network("eo_load_xdiff_R50", cx, z0=50.0).write_touchstone(dir=str(ads_dir))
    make_network("eo_load_ydiff_R50", cy, z0=50.0).write_touchstone(dir=str(ads_dir))

    # =========================================================================
    # Validation Plots
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Z magnitude
    ax = axes[0, 0]
    ax.plot(freq.f/1e6, np.abs(ntwk_2p_50.z[:,0,0]), label='x differential load')
    ax.plot(freq.f/1e6, np.abs(ntwk_2p_50.z[:,1,1]), label='y differential load')
    ax.plot(freq.f/1e6, 1.0 / (2*np.pi*freq.f*c_diff[0,0]), 'k--', alpha=0.5, label='Ideal 1/(wC_x)')
    ax.set_title("Input Impedance Magnitude |Z|")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Ohms")
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True)
    
    # Z phase
    ax = axes[0, 1]
    ax.plot(freq.f/1e6, np.angle(ntwk_2p_50.z[:,0,0], deg=True), label='x differential')
    ax.plot(freq.f/1e6, np.angle(ntwk_2p_50.z[:,1,1], deg=True), label='y differential')
    ax.set_title("Input Impedance Phase")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Degrees")
    ax.set_ylim(-95, -85)
    ax.legend()
    ax.grid(True)
    
    # S11 magnitude (dB)
    ax = axes[1, 0]
    ax.plot(freq.f/1e6, 20*np.log10(np.abs(ntwk_2p_50.s[:,0,0])), label='S11 x (50Ω)')
    ax.plot(freq.f/1e6, 20*np.log10(np.abs(ntwk_2p_50.s[:,1,1])), label='S22 y (50Ω)')
    ax.plot(freq.f/1e6, 20*np.log10(np.abs(ntwk_2p_100.s[:,0,0])), '--', label='S11 x (100Ω)')
    ax.plot(freq.f/1e6, 20*np.log10(np.abs(ntwk_2p_100.s[:,1,1])), '--', label='S22 y (100Ω)')
    ax.set_title("Reflection Magnitude |S11| (dB)")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("dB")
    ax.legend()
    ax.grid(True)
    
    # Current at Vpi
    ax = axes[1, 1]
    Vpi_x = spec["Vpi_x_V"]
    Vpi_y = spec["Vpi_y_V"]
    # I = V / Z
    Ix = Vpi_x / np.abs(ntwk_2p_50.z[:,0,0])
    Iy = Vpi_y / np.abs(ntwk_2p_50.z[:,1,1])
    ax.plot(freq.f/1e6, Ix, label=f'x current @ {Vpi_x:.1f}V')
    ax.plot(freq.f/1e6, Iy, label=f'y current @ {Vpi_y:.1f}V')
    ax.set_title("Peak Current at $V_\\pi$")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Amps")
    ax.legend()
    ax.grid(True)
    
    fig.tight_layout()
    fig.savefig(ads_dir / "validation_plots.png", dpi=150)
    
    # =========================================================================
    # Write README
    # =========================================================================
    readme = """# ADS / SPICE Import Package for Rotating EO Modulator

This directory contains the RF-ready linear load models for the locked 500 µm hybrid geometry.

> **WARNING**: These are strictly small-signal linear load models (capacitive with optional dielectric loss). They do **NOT** model the nonlinear LDMOS PA or the RF matching network itself.

## Touchstone Files (.sNp)
Use these in ADS using the `SnP` component.

### 4-Port Single-Ended (`eo_load_4port_single_ended_*.s4p`)
Reference impedance: 50 ohms.
Port definitions:
1. `+x_full` electrode (left face)
2. `-x_full` electrode (right face)
3. `+y_strip` electrode (top face)
4. `-y_strip` electrode (bottom face)

### 2-Port Differential (`eo_load_2port_diff_*.s2p`)
Reference impedance: 50 ohms or 100 ohms (as marked in filename).
Port definitions:
1. `x` differential pair (driven as $+V/2$ and $-V/2$)
2. `y` differential pair (driven as $+V/2$ and $-V/2$)

Frequency Range: 0.5 MHz to 200 MHz

## SPICE Subcircuits (.sp)
Standard text subcircuits defining pairwise Maxwell capacitances.
`EO_LOAD_4NODE p_x m_x p_y m_y`

For the lossy variants (`tand_1e-4` etc.), parallel resistors are added to simulate the equivalent dielectric conductance evaluated at **100 MHz**. For accurate broadband loss simulation in SPICE, consider using Laplace blocks or replacing the `.sp` with the corresponding Touchstone file in your SPICE simulator.
"""
    with open(ads_dir / "README_ADS_IMPORT.md", "w", encoding="utf-8") as f:
        f.write(readme)
        
    print("Export complete.")

if __name__ == "__main__":
    main()
