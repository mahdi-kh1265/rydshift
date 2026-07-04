"""Hybrid electrode geometry sweep.

Geometry:
  - 3 mm x 3 mm cross-section, 30 mm long z-cut LiNbO3 bar.
  - Existing pair: full-face +x/-x electrodes (entire left/right faces).
  - Added pair: +y/-y strip electrodes with variable setback s from corners.
    Strip width = 3 mm - 2s.
  - Unelectroded top/bottom corner portions are Neumann (insulating).

This script sweeps setback s and evaluates:
  - Full 4x4 Maxwell capacitance matrix C_node
  - 2x2 differential capacitance matrix C_diff (C_xx != C_yy)
  - 2x2 field-per-volt matrix A (A_xx != A_yy)
  - Separate Vpi_x and Vpi_y
  - RF current burden for quadrature drive
  - Field uniformity for each axis
  - Cross-axis leakage
  - Beam offset robustness
  - Field maps
"""
import sys, time
from pathlib import Path

src_path = str(Path(r"c:\Users\khams008\Documents\rydshift\rot_eo_jupyter_handoff\src"))
sys.path.insert(0, src_path)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from skfem import MeshTri
from rot_eo_model.electrostatics import (
    ElectrostaticModel2D,
    gaussian_beam_weights,
    EPSILON_0,
)
from rot_eo_model.tensors import CapacitanceMatrix
from rot_eo_model.eo import half_wave_field_V_per_m
from rot_eo_model.materials import sellmeier_no

out_dir = Path(
    r"C:\Users\khams008\.gemini\antigravity-ide\brain"
    r"\16a39c4f-cd79-46a8-9129-286bf7ecb680"
)

# ── Constants ─────────────────────────────────────────────────────────
CS = 3e-3          # crystal side (m)
HALF = CS / 2
L_EL = 30e-3       # electrode length (m)
EPS_R = 44.0       # in-plane eps for z-cut LiNbO3
W0 = 0.5e-3        # beam waist (m)
LAMBDA0 = 780e-9
N_O = sellmeier_no(LAMBDA0)
R22 = 6.80e-12
EPI = half_wave_field_V_per_m(lambda0_m=LAMBDA0, length_m=L_EL, n_o=N_O, r22_m_per_V=R22)

# Differential transformation: V_node = T @ [Vx, Vy]
# Node order: [+x, -x, +y, -y]
T = 0.5 * np.array([[ 1, 0],
                     [-1, 0],
                     [ 0, 1],
                     [ 0,-1]])

names4 = ['+x', '-x', '+y', '-y']
RES = 150  # mesh resolution per side


# =====================================================================
# Mesh builder: hybrid geometry
# =====================================================================
def build_hybrid_mesh(crystal_size, y_setback, res=150):
    """Build mesh with full-face +x/-x and strip +y/-y electrodes.
    
    +x and -x electrodes cover the entire left/right sides (including corners).
    +y and -y electrodes are centered strips of width (crystal_size - 2*y_setback),
    setback from each corner by y_setback.
    
    Corner nodes belong ONLY to the +x/-x full-face electrodes (they take priority).
    """
    half = crystal_size / 2
    x = np.linspace(-half, half, res)
    y = np.linspace(-half, half, res)
    mesh = MeshTri.init_tensor(x, y)
    tol = crystal_size / (res - 1) * 0.1
    sb = y_setback
    
    # +x/-x: full face including all corners
    # +y/-y: strip, excluding corner regions
    # The strip excludes corners by requiring distance from corner > setback
    mesh = mesh.with_boundaries({
        # Full-face left/right (entire side including corners)
        '+x': lambda p: np.abs(p[0] - half) < tol,
        '-x': lambda p: np.abs(p[0] + half) < tol,
        # Strip top/bottom: must be setback from the left/right edges
        '+y': lambda p: (np.abs(p[1] - half) < tol)
                      & (p[0] > -half + sb + tol) & (p[0] < half - sb - tol),
        '-y': lambda p: (np.abs(p[1] + half) < tol)
                      & (p[0] > -half + sb + tol) & (p[0] < half - sb - tol),
    })
    return mesh


# =====================================================================
# Evaluation grid with optional beam offset
# =====================================================================
def make_eval_grid(offset_x=0.0, offset_y=0.0, n=60):
    """Return (pts, weights) for Gaussian beam at given offset."""
    xg = np.linspace(-CS/2, CS/2, n)
    yg = np.linspace(-CS/2, CS/2, n)
    Xg, Yg = np.meshgrid(xg, yg)
    pts = np.vstack((Xg.ravel(), Yg.ravel()))
    # Beam is centered at (offset_x, offset_y)
    r2 = (pts[0] - offset_x)**2 + (pts[1] - offset_y)**2
    w = np.exp(-2 * r2 / W0**2)
    w /= w.sum()
    return pts, w


# =====================================================================
# Full analysis for one setback + beam position
# =====================================================================
def analyze_config(model, pts, weights, electrode_names=names4):
    """Extract all metrics for a given model and beam position."""
    # ── Capacitance ───────────────────────────────────────────────
    c_prime = model.extract_capacitance_matrix_per_length(electrode_names)
    c_node = c_prime * L_EL
    cap = CapacitanceMatrix(C_node=c_node, T=T)
    c_diff = cap.C_diff
    
    # ── Field-per-volt matrix A ───────────────────────────────────
    # X-drive: Vx=1 => V_node = T @ [1, 0] = [+0.5, -0.5, 0, 0]
    bc_x = {'+x': 0.5, '-x': -0.5, '+y': 0.0, '-y': 0.0}
    ux = model.solve_potential(bc_x)
    E_x = model.extract_electric_field(ux, pts)
    
    # Y-drive: Vy=1 => V_node = T @ [0, 1] = [0, 0, +0.5, -0.5]
    bc_y = {'+x': 0.0, '-x': 0.0, '+y': 0.5, '-y': -0.5}
    uy = model.solve_potential(bc_y)
    E_y = model.extract_electric_field(uy, pts)
    
    # A matrix: [<Ex>, <Ey>] = A @ [Vx, Vy]
    A = np.zeros((2, 2))
    A[0, 0] = np.sum(E_x[0] * weights)  # <Ex> from Vx
    A[1, 0] = np.sum(E_x[1] * weights)  # <Ey> from Vx  (cross-axis)
    A[0, 1] = np.sum(E_y[0] * weights)  # <Ex> from Vy  (cross-axis)
    A[1, 1] = np.sum(E_y[1] * weights)  # <Ey> from Vy
    
    # ── Nonuniformity ─────────────────────────────────────────────
    rms = np.zeros((2, 2))
    rms[0, 0] = np.sqrt(np.sum((E_x[0] - A[0, 0])**2 * weights))
    rms[1, 0] = np.sqrt(np.sum((E_x[1] - A[1, 0])**2 * weights))
    rms[0, 1] = np.sqrt(np.sum((E_y[0] - A[0, 1])**2 * weights))
    rms[1, 1] = np.sqrt(np.sum((E_y[1] - A[1, 1])**2 * weights))
    
    rms_xx_pct = abs(rms[0, 0] / A[0, 0]) * 100 if A[0, 0] != 0 else 0
    rms_yy_pct = abs(rms[1, 1] / A[1, 1]) * 100 if A[1, 1] != 0 else 0
    
    # ── Half-wave voltages ────────────────────────────────────────
    vpi_x = EPI / abs(A[0, 0]) if A[0, 0] != 0 else float('inf')
    vpi_y = EPI / abs(A[1, 1]) if A[1, 1] != 0 else float('inf')
    
    # ── Cross-axis leakage ────────────────────────────────────────
    leak_xy = abs(A[0, 1] / A[0, 0]) if A[0, 0] != 0 else 0  # Ex from Vy / Ex from Vx
    leak_yx = abs(A[1, 0] / A[1, 1]) if A[1, 1] != 0 else 0  # Ey from Vx / Ey from Vy
    
    return {
        'c_node': c_node,
        'c_diff': c_diff,
        'cdiff_xx': c_diff[0, 0],
        'cdiff_yy': c_diff[1, 1],
        'cdiff_xy': c_diff[0, 1],
        'A': A,
        'Axx': A[0, 0], 'Ayy': A[1, 1],
        'Axy': A[0, 1], 'Ayx': A[1, 0],
        'vpi_x': vpi_x, 'vpi_y': vpi_y,
        'rms_xx_pct': rms_xx_pct, 'rms_yy_pct': rms_yy_pct,
        'leak_xy': leak_xy * 100,  # percent
        'leak_yx': leak_yx * 100,
        'ux': ux, 'uy': uy,  # save potential fields for field maps
    }


# =====================================================================
# Main sweep
# =====================================================================
setbacks_um = [100, 250, 500, 750, 1000, 1250]
freqs = [6e6, 10e6, 50e6, 100e6]

pts0, w0 = make_eval_grid(0, 0)

all_results = []

print("=" * 70)
print("HYBRID ELECTRODE SWEEP: full-face +x/-x, strip +y/-y")
print("=" * 70)

for sb_um in setbacks_um:
    sb_m = sb_um * 1e-6
    strip_width_um = 3000 - 2 * sb_um
    
    t0 = time.time()
    mesh = build_hybrid_mesh(CS, sb_m, res=RES)
    model = ElectrostaticModel2D(mesh, lambda p: np.full(p.shape[1], EPS_R))
    m = analyze_config(model, pts0, w0)
    dt = time.time() - t0
    
    # RF currents for each channel at each frequency
    # For pure x-drive at Vpi_x: I_x = omega * C_diff_xx * Vpi_x + omega * C_diff_xy * 0
    # For pure y-drive at Vpi_y: I_y = omega * C_diff_yy * Vpi_y + omega * C_diff_xy * 0
    # For quadrature: V_diff = [Vpi_x, -j*Vpi_y]
    #   I_x(omega) = j*omega*(C_xx*Vpi_x + C_xy*(-j*Vpi_y)) = j*omega*C_xx*Vpi_x + omega*C_xy*Vpi_y
    #   |I_x| = omega * sqrt((C_xx*Vpi_x)^2 + (C_xy*Vpi_y)^2)
    #   I_y(omega) = j*omega*(C_xy*Vpi_x + C_yy*(-j*Vpi_y)) = j*omega*C_xy*Vpi_x + omega*C_yy*Vpi_y
    #   |I_y| = omega * sqrt((C_xy*Vpi_x)^2 + (C_yy*Vpi_y)^2)
    
    rf_currents = {}
    for f in freqs:
        omega = 2 * np.pi * f
        # Single-channel peak currents
        ix_single = omega * abs(m['cdiff_xx']) * m['vpi_x']
        iy_single = omega * abs(m['cdiff_yy']) * m['vpi_y']
        # Quadrature peak currents
        ix_quad = omega * np.sqrt((m['cdiff_xx'] * m['vpi_x'])**2 + (m['cdiff_xy'] * m['vpi_y'])**2)
        iy_quad = omega * np.sqrt((m['cdiff_xy'] * m['vpi_x'])**2 + (m['cdiff_yy'] * m['vpi_y'])**2)
        rf_currents[f] = {
            'ix_single': ix_single, 'iy_single': iy_single,
            'ix_quad': ix_quad, 'iy_quad': iy_quad,
        }
    
    # Beam offset robustness
    offsets = [(0, 0), (100e-6, 0), (200e-6, 0), (0, 100e-6), (0, 200e-6)]
    offset_results = {}
    for ox_m, oy_m in offsets:
        key = f"dx={int(ox_m*1e6)}um,dy={int(oy_m*1e6)}um"
        pts_off, w_off = make_eval_grid(ox_m, oy_m)
        mo = analyze_config(model, pts_off, w_off)
        offset_results[key] = {
            'Axx': mo['Axx'], 'Ayy': mo['Ayy'],
            'vpi_x': mo['vpi_x'], 'vpi_y': mo['vpi_y'],
            'rms_xx_pct': mo['rms_xx_pct'], 'rms_yy_pct': mo['rms_yy_pct'],
            'leak_xy': mo['leak_xy'], 'leak_yx': mo['leak_yx'],
        }
    
    row = {
        'sb_um': sb_um,
        'strip_width_um': strip_width_um,
        **m,
        'rf_currents': rf_currents,
        'offsets': offset_results,
        'model': model,
        'mesh': mesh,
    }
    all_results.append(row)
    
    print(f"\n  Setback: {sb_um} um  Strip width: {strip_width_um} um  ({dt:.1f}s)")
    print(f"    C_diff_xx = {m['cdiff_xx']*1e12:.2f} pF    C_diff_yy = {m['cdiff_yy']*1e12:.2f} pF    C_diff_xy = {m['cdiff_xy']*1e12:.2f} pF")
    print(f"    A_xx = {m['Axx']:.1f}    A_yy = {m['Ayy']:.1f}    A_xy = {m['Axy']:.2f}    A_yx = {m['Ayx']:.2f}")
    print(f"    Vpi_x = {m['vpi_x']:.1f} V    Vpi_y = {m['vpi_y']:.1f} V")
    print(f"    RMS_x = {m['rms_xx_pct']:.2f}%    RMS_y = {m['rms_yy_pct']:.2f}%")
    print(f"    Leak_xy = {m['leak_xy']:.3f}%    Leak_yx = {m['leak_yx']:.3f}%")
    print(f"    I_x@100MHz = {rf_currents[100e6]['ix_single']:.2f} A    I_y@100MHz = {rf_currents[100e6]['iy_single']:.2f} A")


# =====================================================================
# PLOTS
# =====================================================================
print("\nGenerating plots...")

sb_arr = [r['sb_um'] for r in all_results]

# ── Plot 1: C_xx and C_yy vs setback ─────────────────────────────
fig1, ax = plt.subplots(figsize=(8, 5))
cxx = [r['cdiff_xx']*1e12 for r in all_results]
cyy = [r['cdiff_yy']*1e12 for r in all_results]
ax.plot(sb_arr, cxx, 'o-', color='#2196F3', linewidth=2, markersize=7, label='$C_{diff,xx}$ (full-face)')
ax.plot(sb_arr, cyy, 's-', color='#F44336', linewidth=2, markersize=7, label='$C_{diff,yy}$ (strip)')
ax.set_xlabel('Added-pair setback (um)', fontsize=12)
ax.set_ylabel('$C_{diff}$ (pF)', fontsize=12)
ax.set_title('Differential Capacitance: Full-Face vs Strip Electrodes', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig(out_dir / "hybrid_C_vs_setback.png", dpi=150)

# ── Plot 2: Vpi_x and Vpi_y vs setback ───────────────────────────
fig2, ax = plt.subplots(figsize=(8, 5))
vpx = [r['vpi_x'] for r in all_results]
vpy = [r['vpi_y'] for r in all_results]
ax.plot(sb_arr, vpx, 'o-', color='#2196F3', linewidth=2, markersize=7, label='$V_{\\pi,x}$ (full-face drive)')
ax.plot(sb_arr, vpy, 's-', color='#F44336', linewidth=2, markersize=7, label='$V_{\\pi,y}$ (strip drive)')
ax.axhline(700, color='gray', linestyle='--', alpha=0.5, label='700 V limit')
ax.set_xlabel('Added-pair setback (um)', fontsize=12)
ax.set_ylabel('$V_\\pi$ (V)', fontsize=12)
ax.set_title('Half-Wave Voltage: Asymmetry Between Channels', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig(out_dir / "hybrid_Vpi_vs_setback.png", dpi=150)

# ── Plot 3: RF current at 100 MHz ────────────────────────────────
fig3, ax = plt.subplots(figsize=(8, 5))
ix100 = [r['rf_currents'][100e6]['ix_single'] for r in all_results]
iy100 = [r['rf_currents'][100e6]['iy_single'] for r in all_results]
ax.plot(sb_arr, ix100, 'o-', color='#2196F3', linewidth=2, markersize=7, label='$I_{x,peak}$ (full-face)')
ax.plot(sb_arr, iy100, 's-', color='#F44336', linewidth=2, markersize=7, label='$I_{y,peak}$ (strip)')
ax.set_xlabel('Added-pair setback (um)', fontsize=12)
ax.set_ylabel('$I_{peak}$ @ 100 MHz (A)', fontsize=12)
ax.set_title('RF Current Burden per Channel @ 100 MHz', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
fig3.tight_layout()
fig3.savefig(out_dir / "hybrid_I100_vs_setback.png", dpi=150)

# ── Plot 4: Beam offset sensitivity ──────────────────────────────
fig4, axes = plt.subplots(1, 2, figsize=(13, 5))

# Pick 500, 750, 1000 um setbacks for comparison
highlight_sbs = [500, 750, 1000]
colors_hl = {'500': '#4CAF50', '750': '#FF9800', '1000': '#9C27B0'}

for sb_um in highlight_sbs:
    r = [x for x in all_results if x['sb_um'] == sb_um][0]
    off = r['offsets']
    
    # X-offset sensitivity on Vpi_y (the strip channel)
    x_offs = [0, 100, 200]
    vpy_xoff = [off[f"dx={dx}um,dy=0um"]['vpi_y'] for dx in x_offs]
    axes[0].plot(x_offs, vpy_xoff, 'o-', color=colors_hl[str(sb_um)],
                 linewidth=2, markersize=6, label=f's={sb_um} um')
    
    # Y-offset sensitivity on Vpi_y
    y_offs = [0, 100, 200]
    vpy_yoff = [off[f"dx=0um,dy={dy}um"]['vpi_y'] for dy in y_offs]
    axes[1].plot(y_offs, vpy_yoff, 's-', color=colors_hl[str(sb_um)],
                 linewidth=2, markersize=6, label=f's={sb_um} um')

axes[0].set_xlabel('Beam x-offset (um)', fontsize=11)
axes[0].set_ylabel('$V_{\\pi,y}$ (V)', fontsize=11)
axes[0].set_title('$V_{\\pi,y}$ vs x-offset (alignment sensitivity)', fontsize=12)
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].set_xlabel('Beam y-offset (um)', fontsize=11)
axes[1].set_ylabel('$V_{\\pi,y}$ (V)', fontsize=11)
axes[1].set_title('$V_{\\pi,y}$ vs y-offset (alignment sensitivity)', fontsize=12)
axes[1].legend()
axes[1].grid(True, alpha=0.3)
fig4.tight_layout()
fig4.savefig(out_dir / "hybrid_offset_sensitivity.png", dpi=150)

# ── Plot 5: Field maps for 500 um and 750 um setback ─────────────
for sb_target in [500, 750]:
    r = [x for x in all_results if x['sb_um'] == sb_target][0]
    model = r['model']
    mesh = r['mesh']
    
    # Dense field evaluation grid
    Nf = 100
    xf = np.linspace(-CS/2, CS/2, Nf)
    yf = np.linspace(-CS/2, CS/2, Nf)
    Xf, Yf = np.meshgrid(xf, yf)
    pts_f = np.vstack((Xf.ravel(), Yf.ravel()))
    
    # X-drive field
    bc_x = {'+x': 0.5, '-x': -0.5, '+y': 0.0, '-y': 0.0}
    ux = model.solve_potential(bc_x)
    E_xdrive = model.extract_electric_field(ux, pts_f)
    Ex_map = E_xdrive[0].reshape(Nf, Nf)
    
    # Y-drive field
    bc_y = {'+x': 0.0, '-x': 0.0, '+y': 0.5, '-y': -0.5}
    uy = model.solve_potential(bc_y)
    E_ydrive = model.extract_electric_field(uy, pts_f)
    Ey_map = E_ydrive[1].reshape(Nf, Nf)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    
    ext = [-HALF*1e3, HALF*1e3, -HALF*1e3, HALF*1e3]
    
    im1 = ax1.imshow(Ex_map, extent=ext, origin='lower', cmap='RdBu_r', aspect='equal')
    ax1.set_title(f'$E_x$ from x-drive (s={sb_target} um)', fontsize=12)
    ax1.set_xlabel('x (mm)'); ax1.set_ylabel('y (mm)')
    plt.colorbar(im1, ax=ax1, label='V/m per V_diff')
    # Draw beam circle
    theta = np.linspace(0, 2*np.pi, 100)
    ax1.plot(W0*1e3*np.cos(theta), W0*1e3*np.sin(theta), 'k--', linewidth=1, alpha=0.5)
    # Draw electrode boundaries
    sb_mm = sb_target / 1000
    # Full face left/right
    ax1.plot([-1.5, -1.5], [-1.5, 1.5], 'k-', linewidth=3)
    ax1.plot([1.5, 1.5], [-1.5, 1.5], 'k-', linewidth=3)
    # Strip top/bottom
    ax1.plot([-1.5+sb_mm, 1.5-sb_mm], [1.5, 1.5], 'r-', linewidth=3)
    ax1.plot([-1.5+sb_mm, 1.5-sb_mm], [-1.5, -1.5], 'r-', linewidth=3)
    
    im2 = ax2.imshow(Ey_map, extent=ext, origin='lower', cmap='RdBu_r', aspect='equal')
    ax2.set_title(f'$E_y$ from y-drive (s={sb_target} um)', fontsize=12)
    ax2.set_xlabel('x (mm)'); ax2.set_ylabel('y (mm)')
    plt.colorbar(im2, ax=ax2, label='V/m per V_diff')
    ax2.plot(W0*1e3*np.cos(theta), W0*1e3*np.sin(theta), 'k--', linewidth=1, alpha=0.5)
    ax2.plot([-1.5, -1.5], [-1.5, 1.5], 'k-', linewidth=3)
    ax2.plot([1.5, 1.5], [-1.5, 1.5], 'k-', linewidth=3)
    ax2.plot([-1.5+sb_mm, 1.5-sb_mm], [1.5, 1.5], 'r-', linewidth=3)
    ax2.plot([-1.5+sb_mm, 1.5-sb_mm], [-1.5, -1.5], 'r-', linewidth=3)
    
    fig.tight_layout()
    fig.savefig(out_dir / f"hybrid_fieldmap_s{sb_target}.png", dpi=150)

plt.close("all")


# =====================================================================
# MARKDOWN REPORT
# =====================================================================
print("Generating report...")

lines = []
W = lines.append

W("# Hybrid Electrode Geometry Sweep")
W("")
W("> **Geometry**: 3 mm x 3 mm x 30 mm z-cut LiNbO3 bar.")
W("> Full-face electrodes on left/right (+x/-x); strip electrodes with variable setback on top/bottom (+y/-y).")
W("")
W("---")
W("")

# Section 1: Main sweep table
W("## 1. Sweep Results")
W("")
W("| Setback | Strip W | $C_{diff,xx}$ | $C_{diff,yy}$ | $C_{diff,xy}$ | $A_{xx}$ | $A_{yy}$ | $V_{\\pi,x}$ | $V_{\\pi,y}$ | RMS$_x$ | RMS$_y$ | Leak$_{xy}$ | Leak$_{yx}$ |")
W("|---------|---------|--------------|--------------|--------------|---------|---------|------------|------------|--------|--------|------------|------------|")
W("| (um) | (um) | (pF) | (pF) | (pF) | (V/m/V) | (V/m/V) | (V) | (V) | (%) | (%) | (%) | (%) |")

for r in all_results:
    W(f"| {r['sb_um']} | {r['strip_width_um']} | "
      f"{r['cdiff_xx']*1e12:.2f} | {r['cdiff_yy']*1e12:.2f} | {r['cdiff_xy']*1e12:.3f} | "
      f"{r['Axx']:.1f} | {r['Ayy']:.1f} | "
      f"{r['vpi_x']:.1f} | {r['vpi_y']:.1f} | "
      f"{r['rms_xx_pct']:.2f} | {r['rms_yy_pct']:.2f} | "
      f"{r['leak_xy']:.3f} | {r['leak_yx']:.3f} |")

W("")
W("![C_diff vs setback](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/hybrid_C_vs_setback.png)")
W("")
W("![Vpi vs setback](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/hybrid_Vpi_vs_setback.png)")
W("")

# Section 2: C_node matrices for key setbacks
W("---")
W("")
W("## 2. Full $C_{node}$ Matrices (pF)")
W("")
for sb_target in [500, 750, 1000]:
    r = [x for x in all_results if x['sb_um'] == sb_target][0]
    Cn = r['c_node'] * 1e12
    W(f"### Setback = {sb_target} um (strip width = {r['strip_width_um']} um)")
    W("")
    W("| | +x | -x | +y | -y |")
    W("|---|---:|---:|---:|---:|")
    for i, name_i in enumerate(names4):
        vals = " | ".join(f"{Cn[i,j]:.3f}" for j in range(4))
        W(f"| {name_i} | {vals} |")
    W("")

# Section 3: RF Current Budget
W("---")
W("")
W("## 3. RF Current Budget")
W("")
W("### Single-channel peak currents")
W("")
W("| Setback | $I_{x}$ @ 6 | $I_{x}$ @ 10 | $I_{x}$ @ 50 | $I_{x}$ @ 100 | $I_{y}$ @ 6 | $I_{y}$ @ 10 | $I_{y}$ @ 50 | $I_{y}$ @ 100 |")
W("|---------|------------|-------------|-------------|--------------|------------|-------------|-------------|--------------|")
W("| (um) | (A) | (A) | (A) | (A) | (A) | (A) | (A) | (A) |")

for r in all_results:
    rc = r['rf_currents']
    W(f"| {r['sb_um']} | "
      f"{rc[6e6]['ix_single']:.3f} | {rc[10e6]['ix_single']:.3f} | {rc[50e6]['ix_single']:.3f} | {rc[100e6]['ix_single']:.3f} | "
      f"{rc[6e6]['iy_single']:.3f} | {rc[10e6]['iy_single']:.3f} | {rc[50e6]['iy_single']:.3f} | {rc[100e6]['iy_single']:.3f} |")

W("")
W("### Quadrature drive: $V_{diff} = [V_{\\pi,x},\\; -jV_{\\pi,y}]$")
W("")
W("| Setback | $|I_x|$ @ 100 MHz | $|I_y|$ @ 100 MHz |")
W("|---------|-------------------|-------------------|")
for r in all_results:
    rc = r['rf_currents']
    W(f"| {r['sb_um']} | {rc[100e6]['ix_quad']:.3f} | {rc[100e6]['iy_quad']:.3f} |")

W("")
W("![I100 vs setback](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/hybrid_I100_vs_setback.png)")
W("")

# Section 4: Beam offset robustness
W("---")
W("")
W("## 4. Beam Offset Robustness")
W("")

for sb_target in [500, 750, 1000]:
    r = [x for x in all_results if x['sb_um'] == sb_target][0]
    W(f"### Setback = {sb_target} um")
    W("")
    W("| Offset | $A_{xx}$ | $A_{yy}$ | $V_{\\pi,x}$ | $V_{\\pi,y}$ | RMS$_x$ | RMS$_y$ | Leak$_{xy}$ | Leak$_{yx}$ |")
    W("|--------|---------|---------|------------|------------|--------|--------|------------|------------|")
    for key, ov in r['offsets'].items():
        W(f"| {key} | {ov['Axx']:.1f} | {ov['Ayy']:.1f} | {ov['vpi_x']:.1f} | {ov['vpi_y']:.1f} | "
          f"{ov['rms_xx_pct']:.2f} | {ov['rms_yy_pct']:.2f} | {ov['leak_xy']:.3f} | {ov['leak_yx']:.3f} |")
    W("")

W("![Offset sensitivity](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/hybrid_offset_sensitivity.png)")
W("")

# Section 5: Field maps
W("---")
W("")
W("## 5. Field Maps")
W("")
W("Black lines = full-face +x/-x electrodes. Red lines = strip +y/-y electrodes. Dashed circle = 500 um beam waist.")
W("")
for sb_target in [500, 750]:
    W(f"### Setback = {sb_target} um")
    W("")
    W(f"![Field map s={sb_target}](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/hybrid_fieldmap_s{sb_target}.png)")
    W("")

# Section 6: Design recommendation
W("---")
W("")
W("## 6. Design Recommendation")
W("")
W("### Constraint check")
W("")
W("| Setback | $V_{\\pi,x} < 700$ | $V_{\\pi,y} < 700$ | RMS$_x < 7\\%$ | RMS$_y < 7\\%$ | Leak $< 1\\%$ | Gap $\\ge 100$ um |")
W("|---------|:--:|:--:|:--:|:--:|:--:|:--:|")

for r in all_results:
    vx_ok = "PASS" if r['vpi_x'] < 700 else "FAIL"
    vy_ok = "PASS" if r['vpi_y'] < 700 else "FAIL"
    rx_ok = "PASS" if r['rms_xx_pct'] < 7 else "FAIL"
    ry_ok = "PASS" if r['rms_yy_pct'] < 7 else "FAIL"
    lk_ok = "PASS" if max(r['leak_xy'], r['leak_yx']) < 1 else "FAIL"
    gp_ok = "PASS" if r['sb_um'] >= 100 else "FAIL"
    W(f"| {r['sb_um']} | {vx_ok} | {vy_ok} | {rx_ok} | {ry_ok} | {lk_ok} | {gp_ok} |")

W("")

# Find best: lowest total I_peak at 100 MHz while all constraints pass
valid = [r for r in all_results
         if r['vpi_x'] < 700 and r['vpi_y'] < 700
         and r['rms_xx_pct'] < 7 and r['rms_yy_pct'] < 7
         and max(r['leak_xy'], r['leak_yx']) < 1
         and r['sb_um'] >= 100]

if valid:
    # Sort by total RF current at 100 MHz
    best = min(valid, key=lambda r: r['rf_currents'][100e6]['ix_single'] + r['rf_currents'][100e6]['iy_single'])
    rc = best['rf_currents']
    
    W(f"### Recommended: **{best['sb_um']} um setback** (strip width = {best['strip_width_um']} um)")
    W("")
    W("This provides the lowest total RF current while satisfying all constraints:")
    W("")
    W("| Parameter | x-channel (full-face) | y-channel (strip) |")
    W("|-----------|:---------------------:|:-----------------:|")
    W(f"| $C_{{diff}}$ | {best['cdiff_xx']*1e12:.2f} pF | {best['cdiff_yy']*1e12:.2f} pF |")
    W(f"| $A$ | {best['Axx']:.1f} V/m/V | {best['Ayy']:.1f} V/m/V |")
    W(f"| $V_\\pi$ | {best['vpi_x']:.1f} V | {best['vpi_y']:.1f} V |")
    W(f"| RMS nonuniformity | {best['rms_xx_pct']:.2f}% | {best['rms_yy_pct']:.2f}% |")
    W(f"| $I_{{peak}}$ @ 6 MHz | {rc[6e6]['ix_single']:.3f} A | {rc[6e6]['iy_single']:.3f} A |")
    W(f"| $I_{{peak}}$ @ 10 MHz | {rc[10e6]['ix_single']:.3f} A | {rc[10e6]['iy_single']:.3f} A |")
    W(f"| $I_{{peak}}$ @ 50 MHz | {rc[50e6]['ix_single']:.3f} A | {rc[50e6]['iy_single']:.3f} A |")
    W(f"| $I_{{peak}}$ @ 100 MHz | {rc[100e6]['ix_single']:.3f} A | {rc[100e6]['iy_single']:.3f} A |")
    W(f"| Cross-axis leakage | {best['leak_xy']:.3f}% | {best['leak_yx']:.3f}% |")
    
    # Asymmetry note
    vpi_ratio = best['vpi_y'] / best['vpi_x']
    W("")
    W(f"> [!IMPORTANT]")
    W(f"> The y-channel (strip electrode) has $V_{{\\pi,y}}/V_{{\\pi,x}}$ = {vpi_ratio:.2f}.")
    W(f"> This asymmetry must be compensated in the RF drive amplitudes.")
    W(f"> The x-channel capacitance is larger because the full-face electrodes")
    W(f"> include the corner regions that couple to the adjacent (strip) electrodes.")
else:
    W("> [!WARNING]")
    W("> No configuration satisfied all constraints simultaneously.")

report_text = "\n".join(lines) + "\n"
report_path = out_dir / "hybrid_electrode_sweep.md"
with open(report_path, "w", encoding="utf-8") as f:
    f.write(report_text)

print(f"\nReport written to {report_path}")
print("Done.")
