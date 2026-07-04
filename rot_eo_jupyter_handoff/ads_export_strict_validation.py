import json
import numpy as np
import skrf as rf
from pathlib import Path
import re

def main():
    proj_dir = Path(r"c:\Users\khams008\Documents\rydshift\rot_eo_jupyter_handoff")
    ads_dir = proj_dir / "ads_export"
    out_file = Path(r"C:\Users\khams008\.gemini\antigravity-ide\brain\16a39c4f-cd79-46a8-9129-286bf7ecb680\ads_export_strict_validation.md")
    
    lines = []
    W = lines.append
    
    W("# ADS Export Strict File-Level Validation")
    W("")
    
    # 1. First 20 lines of key files
    W("## 1. Touchstone File Headers")
    W("")
    key_files = [
        "eo_load_4port_single_ended_R50.s4p",
        "eo_load_2port_diff_R50.s2p",
        "eo_load_2port_diff_R100.s2p"
    ]
    
    for fname in key_files:
        fpath = ads_dir / fname
        W(f"### `{fname}`")
        W("```")
        try:
            with open(fpath, "r", encoding="utf-8") as f:
                head = [next(f) for _ in range(20)]
            for line in head:
                W(line.rstrip())
        except StopIteration:
            pass # File shorter than 20 lines
        W("```")
        W("")

    # Load JSON source of truth
    with open(proj_dir / "rf_ready_geometry_spec.json", "r") as f:
        spec = json.load(f)
    C_node_src = np.array(spec["C_node_pF"]) * 1e-12
    C_diff_src = np.array(spec["C_diff_pF"]) * 1e-12

    # 2-4. Touchstone Metrics & Reconstruction
    W("## 2. Touchstone Readback & Reconstruction Metrics")
    W("")
    
    def validate_network(fname, C_expected, expected_z0):
        fpath = ads_dir / fname
        ntwk = rf.Network(str(fpath))
        
        ports = ntwk.number_of_ports
        f_min = ntwk.f[0] / 1e6
        f_max = ntwk.f[-1] / 1e6
        n_points = len(ntwk.f)
        z0_read = ntwk.z0[0,0] if ntwk.z0.ndim > 1 else ntwk.z0[0]
        
        # Reciprocity
        # max_err = max |S_ij - S_ji|
        S = ntwk.s
        recip_err = np.max(np.abs(S - np.transpose(S, axes=(0, 2, 1))))
        
        # Passivity
        # Max singular value over all frequencies
        s_svd = np.linalg.svd(S, compute_uv=False)
        max_sv = np.max(s_svd)
        is_passive = max_sv <= 1.0 + 1e-6
        
        # Recover Y from S
        # Y = ntwk.y is already computed by skrf
        Y = ntwk.y
        omega = 2 * np.pi * ntwk.f[:, np.newaxis, np.newaxis]
        
        # C_recovered = imag(Y) / omega
        # Guard against zero frequency, but Touchstone starts at 0.5 MHz here
        C_rec = np.imag(Y) / omega
        
        # Compare to expected
        # Max relative error
        abs_err = np.abs(C_rec - C_expected[np.newaxis, :, :])
        # Mask zeros in C_expected to avoid division by zero
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_err = abs_err / np.abs(C_expected[np.newaxis, :, :])
            rel_err = np.where(np.abs(C_expected)[np.newaxis, :, :] < 1e-15, 0.0, rel_err)
        
        max_rel_err = np.max(rel_err)
        
        W(f"### `{fname}`")
        W(f"- **Ports**: {ports}")
        W(f"- **Frequency**: {f_min:.1f} to {f_max:.1f} MHz ({n_points} pts)")
        W(f"- **$Z_0$ Read**: {z0_read} $\\Omega$")
        W(f"- **Reciprocity Error**: {recip_err:.3e}")
        W(f"- **Max Singular Value**: {max_sv:.6f} (Passive: {is_passive})")
        W(f"- **Max Capacitance Reconstruction Error**: {max_rel_err*100:.4e} %")
        W("")

    validate_network("eo_load_4port_single_ended_R50.s4p", C_node_src, 50.0)
    validate_network("eo_load_2port_diff_R50.s2p", C_diff_src, 50.0)
    validate_network("eo_load_2port_diff_R100.s2p", C_diff_src, 100.0)
    validate_network("eo_load_4port_single_ended_R50_tand_0.001.s4p", C_node_src, 50.0)
    validate_network("eo_load_2port_diff_R50_tand_0.001.s2p", C_diff_src, 50.0)

    # 5. SPICE Subcircuit Validation
    W("## 5. SPICE Subcircuit Validation")
    W("")
    W("Extracting capacitors from `eo_load_4node.sp` and reconstructing $C_{node}$:")
    W("")
    
    sp_path = ads_dir / "eo_load_4node.sp"
    with open(sp_path, "r") as f:
        sp_content = f.read()
    
    node_map = {"p_x": 0, "m_x": 1, "p_y": 2, "m_y": 3, "0": None}
    
    # Parse capacitors
    caps = []
    C_recon = np.zeros((4, 4))
    
    W("```spice")
    for line in sp_content.splitlines():
        line = line.strip()
        if line.startswith("C"):
            # e.g., C1 p_x m_x 2.943527e-12
            W(line)
            parts = line.split()
            n1, n2, val_str = parts[1], parts[2], parts[3]
            val = float(val_str)
            i1 = node_map.get(n1)
            i2 = node_map.get(n2)
            
            if i1 is not None and i2 is not None:
                # Pairwise capacitor C_ij = -C_node[i,j]
                C_recon[i1, i2] = -val
                C_recon[i2, i1] = -val
                # Diagonal contribution: node self-capacitance must include all connected caps
                C_recon[i1, i1] += val
                C_recon[i2, i2] += val
            elif i1 is not None and i2 is None:
                # Capacitor to ground
                C_recon[i1, i1] += val
            elif i2 is not None and i1 is None:
                C_recon[i2, i2] += val
    W("```")
    W("")
    
    abs_err_spice = np.abs(C_recon - C_node_src)
    max_abs_spice = np.max(abs_err_spice)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_err_spice = abs_err_spice / np.abs(C_node_src)
        rel_err_spice[np.abs(C_node_src) < 1e-15] = 0.0
    max_rel_spice = np.max(rel_err_spice)
    
    W(f"- **Max Absolute Error**: {max_abs_spice:.3e} F")
    W(f"- **Max Relative Error**: {max_rel_spice*100:.4e} %")
    W("")
    W(f"SPICE construction matches source matrix: **{max_abs_spice < 1e-15}**")
    W("")

    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        
if __name__ == "__main__":
    main()
