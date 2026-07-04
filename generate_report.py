"""Numerical audit report generator — full-face electrode geometry.

Crystal: z-cut LiNbO3, 3 mm x 3 mm x 30 mm bar.
Electrodes: full-face on all four sides (each 3 mm x 30 mm).
Optical propagation: along the 30 mm z-axis.
"""
import sys
from pathlib import Path

src_path = str(Path(r"c:\Users\khams008\Documents\rydshift\rot_eo_jupyter_handoff\src"))
sys.path.insert(0, src_path)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rot_eo_model.constants import DEFAULTS
from rot_eo_model.materials import sellmeier_no, sellmeier_ne, linbo3_dielectric_tensor
from rot_eo_model.eo import (
    half_wave_field_V_per_m, vpi_from_gap_V,
    transverse_eta_matrix, eigen_indices, axis_angle_rad,
)
from rot_eo_model.electrostatics import (
    build_parallel_plate_mesh,
    build_full_face_electrode_mesh,
    ElectrostaticModel2D,
    extract_field_per_volt_matrix,
    EPSILON_0,
)
from rot_eo_model.tensors import CapacitanceMatrix

out_dir = Path(
    r"C:\Users\khams008\.gemini\antigravity-ide\brain"
    r"\16a39c4f-cd79-46a8-9129-286bf7ecb680"
)

# ── Physical constants ────────────────────────────────────────────────
crystal_size = 3e-3          # 3 mm square cross-section
L_electrode  = 30e-3         # electrode length along z
gap          = crystal_size  # electrode-to-electrode gap = crystal side
w0           = 0.5e-3        # Gaussian beam waist radius (500 um)
lambda0      = 780e-9

# ── 1. Material properties ────────────────────────────────────────────
n_o  = sellmeier_no(lambda0)
n_e  = sellmeier_ne(lambda0)
r22  = DEFAULTS.r22_m_per_V
dt   = linbo3_dielectric_tensor(10e6).tensor
eps_11 = dt[0, 0]
eps_33 = dt[2, 2]

# ── 2. EO sensitivity (ideal parallel plate) ─────────────────────────
epi_30 = half_wave_field_V_per_m(lambda0_m=lambda0, length_m=30e-3, n_o=n_o, r22_m_per_V=r22)
vpi_30 = vpi_from_gap_V(gap, lambda0_m=lambda0, length_m=30e-3, n_o=n_o, r22_m_per_V=r22)
epi_3  = half_wave_field_V_per_m(lambda0_m=lambda0, length_m=3e-3,  n_o=n_o, r22_m_per_V=r22)
vpi_3  = vpi_from_gap_V(gap, lambda0_m=lambda0, length_m=3e-3,  n_o=n_o, r22_m_per_V=r22)

# ── 3. Parallel-plate FEM benchmark ──────────────────────────────────
pp_width = 10e-3
pp_gap   = 1e-3
mesh_pp  = build_parallel_plate_mesh(pp_width, pp_gap, res=30)
model_pp = ElectrostaticModel2D(mesh_pp, lambda p: np.full(p.shape[1], eps_11))
c_pp_fem = model_pp.extract_capacitance_matrix_per_length(["top", "bottom"])[0, 0]
c_pp_ana = EPSILON_0 * eps_11 * pp_width / pp_gap
pct_err  = abs(c_pp_fem - c_pp_ana) / c_pp_ana * 100

# ── 4. Four-electrode model — FULL-FACE ──────────────────────────────
mesh_4el  = build_full_face_electrode_mesh(crystal_size, res=60)
model_4el = ElectrostaticModel2D(mesh_4el, lambda p: np.full(p.shape[1], eps_11))
electrode_names = ["+x", "-x", "+y", "-y"]

c_node_prime = model_4el.extract_capacitance_matrix_per_length(electrode_names)
c_node       = c_node_prime * L_electrode

T = 0.5 * np.array([
    [ 1,  0],
    [-1,  0],
    [ 0,  1],
    [ 0, -1],
])
cap = CapacitanceMatrix(C_node=c_node, T=T)
c_diff = cap.C_diff
cross_coupling = c_diff[0, 1] / np.sqrt(c_diff[0, 0] * c_diff[1, 1])

# ── 5. Field-per-volt matrix A ───────────────────────────────────────
Ng = 80
xg = np.linspace(-crystal_size / 2, crystal_size / 2, Ng)
yg = np.linspace(-crystal_size / 2, crystal_size / 2, Ng)
Xg, Yg = np.meshgrid(xg, yg)
points  = np.vstack((Xg.ravel(), Yg.ravel()))

field_metrics = extract_field_per_volt_matrix(model_4el, electrode_names, points, w0)
A   = field_metrics.A
rms = field_metrics.rms_nonuniformity

# ── 6. Cross-axis leakage ────────────────────────────────────────────
leakage_xy = abs(A[0, 1] / A[0, 0]) if A[0, 0] != 0 else 0
leakage_yx = abs(A[1, 0] / A[1, 1]) if A[1, 1] != 0 else 0

# ── 7. FEM-based Vpi ─────────────────────────────────────────────────
# Vpi_ideal  = Epi * gap   (parallel-plate assumption)
# Vpi_FEM    = Epi / |A_xx|  (accounts for actual field distribution)
vpi_ideal = epi_30 * gap
vpi_fem   = epi_30 / abs(A[0, 0])

# ── 8. RF current budget using Vpi_FEM ────────────────────────────────
freqs = [6e6, 10e6, 50e6, 100e6]
c_diff_xx = c_diff[0, 0]
currents_fem = [2 * np.pi * f * c_diff_xx * vpi_fem for f in freqs]

# ── 9. Plots ──────────────────────────────────────────────────────────
# 9a. EO rotation sweep
phi_drive = np.linspace(0, 2 * np.pi, 200)
retardances, angles = [], []
for ex, ey in zip(epi_30 * np.cos(phi_drive), epi_30 * np.sin(phi_drive)):
    eta = transverse_eta_matrix(ex, ey, 0.0, n_o=n_o, pockels=DEFAULTS.pockels)
    n_eig, vecs = eigen_indices(eta)
    retardances.append(2 * np.pi / lambda0 * 30e-3 * abs(n_eig[1] - n_eig[0]))
    angles.append(axis_angle_rad(vecs))

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(np.degrees(phi_drive), np.array(retardances) / np.pi)
ax1.set_xlabel("Drive Phase (deg)")
ax1.set_ylabel(r"Retardance / $\pi$")
ax1.set_title("Retardance Ripple (exact tensor)")
ax1.grid(True)
ax2.plot(np.degrees(phi_drive), np.degrees(angles))
ax2.set_xlabel("Drive Phase (deg)")
ax2.set_ylabel(r"Eigenaxis Angle $\theta$ (deg)")
ax2.set_title("Eigenaxis Rotation")
ax2.grid(True)
plt.tight_layout()
fig1.savefig(out_dir / "eo_rotation.png", dpi=150)

# 9b. X-drive potential & field map
bc_x = {"+x": 0.5, "-x": -0.5, "+y": 0.0, "-y": 0.0}
ux = model_4el.solve_potential(bc_x)
Ex_vec = model_4el.extract_electric_field(ux, points)

fig2, ax = plt.subplots(figsize=(6, 5))
c = ax.tricontourf(mesh_4el.p[0] * 1e3, mesh_4el.p[1] * 1e3, ux, levels=30, cmap="RdBu_r")
fig2.colorbar(c, ax=ax, label="Potential (V) for 1 V diff x-drive")
step = 5
ax.quiver(
    points[0, ::step] * 1e3, points[1, ::step] * 1e3,
    Ex_vec[0, ::step], Ex_vec[1, ::step],
    color="k", alpha=0.4, scale=5000,
)
ax.set_title("Full-face X-drive: Potential & E-field")
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_aspect("equal")
plt.tight_layout()
fig2.savefig(out_dir / "field_map.png", dpi=150)

# 9c. Y-drive potential & field map
bc_y = {"+x": 0.0, "-x": 0.0, "+y": 0.5, "-y": -0.5}
uy = model_4el.solve_potential(bc_y)
Ey_vec = model_4el.extract_electric_field(uy, points)

fig3, ax3 = plt.subplots(figsize=(6, 5))
c3 = ax3.tricontourf(mesh_4el.p[0] * 1e3, mesh_4el.p[1] * 1e3, uy, levels=30, cmap="RdBu_r")
fig3.colorbar(c3, ax=ax3, label="Potential (V) for 1 V diff y-drive")
ax3.quiver(
    points[0, ::step] * 1e3, points[1, ::step] * 1e3,
    Ey_vec[0, ::step], Ey_vec[1, ::step],
    color="k", alpha=0.4, scale=5000,
)
ax3.set_title("Full-face Y-drive: Potential & E-field")
ax3.set_xlabel("x (mm)")
ax3.set_ylabel("y (mm)")
ax3.set_aspect("equal")
plt.tight_layout()
fig3.savefig(out_dir / "field_map_y.png", dpi=150)

plt.close("all")

# ── Report ────────────────────────────────────────────────────────────
report = f"""# Numerical Audit Report — Full-Face Electrode Geometry

## Device Under Test
- **Crystal**: z-cut congruent LiNbO\\u2083, 3 mm x 3 mm x 30 mm bar
- **Optical propagation**: along the 30 mm z-axis
- **Electrodes**: full-face on all four sides (each 3 mm wide x 30 mm long)
- **Electrode gap** (face-to-face): 3 mm (= crystal side length)

---

## 1. Material Properties

| Parameter | Value |
|-----------|-------|
| $n_o$ (780 nm) | {n_o:.4f} |
| $n_e$ (780 nm) | {n_e:.4f} |
| $r_{{22}}$ | {r22:.2e} m/V |
| $\\varepsilon_{{11}}$ (RF, 10 MHz) | {eps_11:.1f} |
| $\\varepsilon_{{33}}$ (RF, 10 MHz) | {eps_33:.1f} |

---

## 2. Electro-Optic Sensitivity

| Quantity | L = 30 mm | L = 3 mm |
|----------|-----------|----------|
| $E_\\pi$ (V/cm) | {epi_30/100:.1f} | {epi_3/100:.1f} |
| $V_{{\\pi,ideal}}$ = $E_\\pi \\cdot d$ (V) | {vpi_30:.1f} | {vpi_3:.1f} |

---

## 3. Parallel-Plate FEM Benchmark

| Quantity | Value |
|----------|-------|
| Geometry | {pp_width*1e3:.0f} mm x {pp_gap*1e3:.0f} mm, $\\varepsilon_r$ = {eps_11:.0f} |
| Analytical $C'$ | {c_pp_ana*1e12:.3f} pF/m |
| FEM $C'$ | {c_pp_fem*1e12:.3f} pF/m |
| Error | {pct_err:.2f}% |

---

## 4. Four-Electrode Maxwell Capacitance (Full-Face)

### Boundary conditions used
For each column $j$ of $C'_{{node}}$: electrode $j$ is set to 1 V, all other electrodes to 0 V.
Row $i$ reports the charge per unit length $Q'_i$ on electrode $i$.
Electrode ordering: [+x, \\u2212x, +y, \\u2212y].

### $C'_{{node}}$ (pF/m)
```
{np.array2string(c_node_prime * 1e12, precision=2, suppress_small=True)}
```

### $C_{{node}}$ total for $L_{{electrode}}$ = 30 mm (pF)
```
{np.array2string(c_node * 1e12, precision=2, suppress_small=True)}
```

### $C_{{diff}} = T^T C_{{node}} T$ (pF)
```
{np.array2string(c_diff * 1e12, precision=3, suppress_small=True)}
```

- **$C_{{diff,xx}}$**: {c_diff[0,0]*1e12:.3f} pF
- **$C_{{diff,yy}}$**: {c_diff[1,1]*1e12:.3f} pF
- **Cross-coupling** $C_{{xy}} / \\sqrt{{C_{{xx}} C_{{yy}}}}$: {cross_coupling:.3e}

---

## 5. Field-per-Volt Matrix $A$

Gaussian-beam-weighted average field (beam waist $w_0$ = {w0*1e6:.0f} $\\mu$m) for 1 V differential drive.

### Boundary conditions
- **X-drive**: $V_{{+x}}$ = +0.5 V, $V_{{-x}}$ = \\u22120.5 V, $V_{{+y}}$ = $V_{{-y}}$ = 0 V
- **Y-drive**: $V_{{+y}}$ = +0.5 V, $V_{{-y}}$ = \\u22120.5 V, $V_{{+x}}$ = $V_{{-x}}$ = 0 V

### $A$ matrix (V/m per V_diff)
```
A = {np.array2string(A, precision=2, suppress_small=True)}
```

- $A_{{xx}}$ = {A[0,0]:.2f} V/m per V
- $A_{{yy}}$ = {A[1,1]:.2f} V/m per V

---

## 6. Cross-Axis Field Leakage

| Ratio | Value |
|-------|-------|
| $|A_{{xy}} / A_{{xx}}|$ | {leakage_xy:.3e} |
| $|A_{{yx}} / A_{{yy}}|$ | {leakage_yx:.3e} |

---

## 7. Gaussian Beam Field Nonuniformity

RMS deviation from the beam-weighted mean field, normalized:

| Component | RMS (V/m per V) | Relative to mean |
|-----------|-----------------|------------------|
| $E_x$ from x-drive | {rms[0,0]:.2f} | {abs(rms[0,0]/A[0,0])*100:.1f}% |
| $E_y$ from y-drive | {rms[1,1]:.2f} | {abs(rms[1,1]/A[1,1])*100:.1f}% |

---

## 8. Vpi Comparison: Ideal vs FEM

| Quantity | Value |
|----------|-------|
| $V_{{\\pi,ideal}}$ = $E_\\pi \\cdot d$ | {vpi_ideal:.1f} V |
| $V_{{\\pi,FEM}}$  = $E_\\pi / |A_{{xx}}|$ | {vpi_fem:.1f} V |
| Ratio $V_{{\\pi,FEM}} / V_{{\\pi,ideal}}$ | {vpi_fem/vpi_ideal:.3f} |

---

## 9. Required RF Current (using $V_{{\\pi,FEM}}$ = {vpi_fem:.1f} V)

$I_{{peak}} = 2\\pi f \\cdot C_{{diff,xx}} \\cdot V_{{\\pi,FEM}}$

| Frequency | $I_{{peak}}$ (A) |
|-----------|-------------------|
| 6 MHz | {currents_fem[0]:.3f} |
| 10 MHz | {currents_fem[1]:.3f} |
| 50 MHz | {currents_fem[2]:.3f} |
| 100 MHz | {currents_fem[3]:.3f} |

---

## 10. Plots

![EO Rotation](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/eo_rotation.png)

![X-drive Field Map](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/field_map.png)

![Y-drive Field Map](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/field_map_y.png)
"""

report_path = out_dir / "numerical_audit_report.md"
with open(report_path, "w", encoding="utf-8") as f:
    f.write(report)

print(f"Report written to {report_path}")
print(f"  C'_node diag: {c_node_prime[0,0]*1e12:.2f} pF/m")
print(f"  C_diff_xx:    {c_diff[0,0]*1e12:.3f} pF")
print(f"  A_xx:         {A[0,0]:.2f} V/m per V")
print(f"  Vpi_ideal:    {vpi_ideal:.1f} V")
print(f"  Vpi_FEM:      {vpi_fem:.1f} V")
