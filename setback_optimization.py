"""Electrode setback optimization sweep.

Sweeps electrode setback to find the optimal trade-off between capacitance
reduction, Vpi penalty, and field uniformity.
"""
import sys
from pathlib import Path

src_path = str(Path(r"c:\Users\khams008\Documents\rydshift\rot_eo_jupyter_handoff\src"))
sys.path.insert(0, src_path)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from skfem import MeshTri
from rot_eo_model.electrostatics import (
    ElectrostaticModel2D,
    gaussian_beam_weights,
)
from rot_eo_model.tensors import CapacitanceMatrix
from rot_eo_model.eo import half_wave_field_V_per_m
from rot_eo_model.materials import sellmeier_no

out_dir = Path(
    r"C:\Users\khams008\.gemini\antigravity-ide\brain"
    r"\16a39c4f-cd79-46a8-9129-286bf7ecb680"
)

# Constants
CS = 3e-3
L_EL = 30e-3
EPS_R = 44.0
W0 = 0.5e-3
LAMBDA0 = 780e-9
N_O = sellmeier_no(LAMBDA0)
R22 = 6.80e-12
EPI = half_wave_field_V_per_m(lambda0_m=LAMBDA0, length_m=L_EL, n_o=N_O, r22_m_per_V=R22)

T = 0.5 * np.array([[ 1, 0], [-1, 0], [ 0, 1], [ 0,-1]])
names4 = ['+x', '-x', '+y', '-y']

Ng = 60
xg = np.linspace(-CS/2, CS/2, Ng)
yg = np.linspace(-CS/2, CS/2, Ng)
Xg, Yg = np.meshgrid(xg, yg)
pts = np.vstack((Xg.ravel(), Yg.ravel()))
weights = gaussian_beam_weights(pts, W0)

def build_setback_mesh(crystal_size, setback, res=100):
    half = crystal_size / 2
    x = np.linspace(-half, half, res)
    y = np.linspace(-half, half, res)
    mesh = MeshTri.init_tensor(x, y)
    tol = crystal_size / (res - 1) * 0.1
    sb = setback
    mesh = mesh.with_boundaries({
        '+x': lambda p: (np.abs(p[0] - half) < tol)
                      & (p[1] > -half + sb + tol) & (p[1] < half - sb - tol),
        '-x': lambda p: (np.abs(p[0] + half) < tol)
                      & (p[1] > -half + sb + tol) & (p[1] < half - sb - tol),
        '+y': lambda p: (np.abs(p[1] - half) < tol)
                      & (p[0] > -half + sb + tol) & (p[0] < half - sb - tol),
        '-y': lambda p: (np.abs(p[1] + half) < tol)
                      & (p[0] > -half + sb + tol) & (p[0] < half - sb - tol),
    })
    return mesh

def get_max_field(model, ux):
    # E-field at all mesh nodes
    E_nodes = model.extract_electric_field(ux, model.mesh.p)
    E_norm = np.sqrt(E_nodes[0]**2 + E_nodes[1]**2)
    return np.max(E_norm)

def optimize_setback():
    setbacks_um = [100, 150, 200, 250, 350, 500, 650, 750, 900, 1000, 1250]
    RES = 150  # good balance of resolution and speed
    
    results = []
    
    print(f"Sweeping setback from {setbacks_um[0]} to {setbacks_um[-1]} um...")
    
    for sb_um in setbacks_um:
        sb_m = sb_um * 1e-6
        el_width_um = 3000 - 2 * sb_um
        
        mesh = build_setback_mesh(CS, sb_m, res=RES)
        model = ElectrostaticModel2D(mesh, lambda p: np.full(p.shape[1], EPS_R))
        
        # Capacitance
        c_prime = model.extract_capacitance_matrix_per_length(names4)
        c_node = c_prime * L_EL
        cap = CapacitanceMatrix(C_node=c_node, T=T)
        c_diff = cap.C_diff
        cdiff_xx = c_diff[0, 0]
        
        # Field evaluation for X-drive
        bc_x = {'+x': 0.5, '-x': -0.5, '+y': 0.0, '-y': 0.0}
        ux = model.solve_potential(bc_x)
        E_vec = model.extract_electric_field(ux, pts)
        
        Axx = float(np.sum(E_vec[0] * weights))
        Ayx = float(np.sum(E_vec[1] * weights)) # cross-axis leakage (Ey from Vx)
        
        rms_xx = float(np.sqrt(np.sum((E_vec[0] - Axx)**2 * weights)))
        rms_rel = abs(rms_xx / Axx) * 100 if Axx != 0 else 0
        leakage_pct = abs(Ayx / Axx) * 100 if Axx != 0 else 0
        
        # Max field enhancement
        max_E = get_max_field(model, ux)
        field_enhancement = max_E / abs(Axx) if Axx != 0 else 0
        
        vpi_fem = EPI / abs(Axx) if Axx != 0 else float('inf')
        cv_product = cdiff_xx * vpi_fem
        ipeak_100 = 2 * np.pi * 100e6 * cdiff_xx * vpi_fem
        
        res_dict = {
            'sb_um': sb_um,
            'el_width_um': el_width_um,
            'cdiff_pF': cdiff_xx * 1e12,
            'Axx': Axx,
            'vpi_fem': vpi_fem,
            'cv_nC': cv_product * 1e9,
            'ipeak_100': ipeak_100,
            'rms_pct': rms_rel,
            'leakage_pct': leakage_pct,
            'max_E_enh': field_enhancement
        }
        results.append(res_dict)
        print(f"  sb={sb_um:4d}um  W={el_width_um:4d}um  C_diff={res_dict['cdiff_pF']:5.2f}pF  "
              f"Vpi={vpi_fem:6.1f}V  I100={ipeak_100:5.2f}A  "
              f"RMS={rms_rel:4.1f}%  Leak={leakage_pct:4.2f}%")
              
    return results

if __name__ == "__main__":
    results = optimize_setback()
    
    # ---------------------------------------------------------
    # Plots
    # ---------------------------------------------------------
    sb_arr = [r['sb_um'] for r in results]
    cdiff_arr = [r['cdiff_pF'] for r in results]
    vpi_arr = [r['vpi_fem'] for r in results]
    i100_arr = [r['ipeak_100'] for r in results]
    
    fig, ax1 = plt.subplots(figsize=(8, 5))
    
    color1 = '#2196F3'
    ax1.set_xlabel('Corner Setback (um)')
    ax1.set_ylabel('$C_{diff,xx}$ (pF)', color=color1)
    ln1 = ax1.plot(sb_arr, cdiff_arr, 'o-', color=color1, label='$C_{diff}$ (pF)')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, alpha=0.3)
    
    ax2 = ax1.twinx()
    color2 = '#FF9800'
    ax2.set_ylabel('$V_{\\pi}$ (V)', color=color2)
    ln2 = ax2.plot(sb_arr, vpi_arr, 's-', color=color2, label='$V_{\\pi}$ (V)')
    ax2.tick_params(axis='y', labelcolor=color2)
    
    # Optional 3rd axis for current (or can just plot on same twinx if scaled)
    # Let's add current text labels or a separate plot, or normalize
    # To keep it simple, we plot Cdiff and Vpi on dual axes, and create a separate subplot for current
    fig.tight_layout()
    fig.savefig(out_dir / "opt_tradeoff_C_Vpi.png", dpi=150)
    
    # Plot 2: Current
    fig2, ax = plt.subplots(figsize=(8, 4))
    ax.plot(sb_arr, i100_arr, '^-', color='#9C27B0', markersize=8)
    ax.set_xlabel('Corner Setback (um)')
    ax.set_ylabel('$I_{peak}$ @ 100 MHz (A)')
    ax.set_title('RF Current vs Setback')
    ax.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(out_dir / "opt_current.png", dpi=150)
    
    plt.close('all')
    
    # ---------------------------------------------------------
    # Best Setback Selection
    # ---------------------------------------------------------
    # Constraints:
    # - Vpi_FEM < 700 V
    # - field nonuniformity < 7%
    # - cross-axis leakage < 1%
    # - fabrication gap >= 100 µm
    
    valid_configs = []
    for r in results:
        if r['vpi_fem'] < 700 and r['rms_pct'] < 7.0 and r['leakage_pct'] < 1.0 and r['sb_um'] >= 100:
            valid_configs.append(r)
            
    if valid_configs:
        # Best is typically lowest RF current / C*Vpi product
        best = min(valid_configs, key=lambda x: x['ipeak_100'])
    else:
        best = None
        
    # ---------------------------------------------------------
    # Markdown Report
    # ---------------------------------------------------------
    lines = []
    W = lines.append
    
    W("# Electrode Setback Optimization Report")
    W("")
    W("Sweeping electrode setback from 100 µm to 1250 µm to minimize RF current while maintaining field uniformity.")
    W("")
    W("## Sweep Results")
    W("")
    W("| Setback | El. Width | $C_{diff,xx}$ | $A_{xx}$ | $V_\\pi$ | $C \\cdot V_\\pi$ | $I_{peak}$ @ 100MHz | Nonuniformity | Leakage | Max Field Enh. |")
    W("|---------|-----------|---------------|----------|--------|-------------------|-------------------|---------------|---------|----------------|")
    W("| (µm) | (µm) | (pF) | (V/m/V) | (V) | (nC) | (A) | (RMS %) | (%) | (Max/Avg) |")
    for r in results:
        W(f"| {r['sb_um']} | {r['el_width_um']} | {r['cdiff_pF']:.2f} | {r['Axx']:.1f} | {r['vpi_fem']:.1f} | {r['cv_nC']:.2f} | {r['ipeak_100']:.2f} | {r['rms_pct']:.2f} | {r['leakage_pct']:.3f} | {r['max_E_enh']:.1f} |")
        
    W("")
    W("![Tradeoff](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/opt_tradeoff_C_Vpi.png)")
    W("![Current](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/opt_current.png)")
    W("")
    W("## Constraint Check & Optimal Selection")
    W("")
    W("**Constraints:**")
    W("- $V_\\pi < 700$ V")
    W("- Nonuniformity < 7%")
    W("- Cross-axis leakage < 1%")
    W("- Setback $\\ge$ 100 µm")
    W("")
    
    if best:
        W("### Selected Best Configuration")
        W(f"**Optimal Setback: {best['sb_um']} µm** (Electrode width: {best['el_width_um']} µm)")
        W("")
        W("This provides the lowest RF current while satisfying all constraints:")
        W(f"- **$C_{{diff,xx}}$:** {best['cdiff_pF']:.2f} pF")
        W(f"- **$V_\\pi$:** {best['vpi_fem']:.1f} V")
        W(f"- **$I_{{peak}}$ @ 100 MHz:** {best['ipeak_100']:.2f} A")
        W(f"- **Nonuniformity:** {best['rms_pct']:.2f}% (Limit: 7%)")
        W(f"- **Cross-axis Leakage:** {best['leakage_pct']:.3f}% (Limit: 1%)")
        W(f"- **Field Enhancement:** {best['max_E_enh']:.1f}$\\times$ (Edge vs Center)")
    else:
        W("No configuration satisfied all constraints.")
        
    with open(out_dir / "setback_optimization.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        
    print("Report written.")
