"""Decision report generator: 500 um vs 750 um added-electrode setback."""
import sys, time
from pathlib import Path

src_path = str(Path(r"c:\Users\khams008\Documents\rydshift\rot_eo_jupyter_handoff\src"))
sys.path.insert(0, src_path)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from skfem import MeshTri
from rot_eo_model.electrostatics import ElectrostaticModel2D, gaussian_beam_weights, EPSILON_0
from rot_eo_model.tensors import CapacitanceMatrix
from rot_eo_model.eo import half_wave_field_V_per_m
from rot_eo_model.materials import sellmeier_no

out_dir = Path(
    r"C:\Users\khams008\.gemini\antigravity-ide\brain"
    r"\16a39c4f-cd79-46a8-9129-286bf7ecb680"
)

# Constants
CS = 3e-3
HALF = CS / 2
L_EL = 30e-3
EPS_R = 44.0
W0 = 0.5e-3
LAMBDA0 = 780e-9
N_O = sellmeier_no(LAMBDA0)
R22 = 6.80e-12
EPI = half_wave_field_V_per_m(lambda0_m=LAMBDA0, length_m=L_EL, n_o=N_O, r22_m_per_V=R22)

T = 0.5 * np.array([[ 1, 0], [-1, 0], [ 0, 1], [ 0,-1]])
names4 = ['+x', '-x', '+y', '-y']
RES = 150

def build_hybrid_mesh(crystal_size, y_setback, res=150):
    half = crystal_size / 2
    x = np.linspace(-half, half, res)
    y = np.linspace(-half, half, res)
    mesh = MeshTri.init_tensor(x, y)
    tol = crystal_size / (res - 1) * 0.1
    sb = y_setback
    mesh = mesh.with_boundaries({
        '+x': lambda p: np.abs(p[0] - half) < tol,
        '-x': lambda p: np.abs(p[0] + half) < tol,
        '+y': lambda p: (np.abs(p[1] - half) < tol) & (p[0] > -half + sb + tol) & (p[0] < half - sb - tol),
        '-y': lambda p: (np.abs(p[1] + half) < tol) & (p[0] > -half + sb + tol) & (p[0] < half - sb - tol),
    })
    return mesh

def make_eval_grid(offset_x=0.0, offset_y=0.0, n=60):
    xg = np.linspace(-CS/2, CS/2, n)
    yg = np.linspace(-CS/2, CS/2, n)
    Xg, Yg = np.meshgrid(xg, yg)
    pts = np.vstack((Xg.ravel(), Yg.ravel()))
    r2 = (pts[0] - offset_x)**2 + (pts[1] - offset_y)**2
    w = np.exp(-2 * r2 / W0**2)
    w /= w.sum()
    return pts, w

def analyze_config(model, pts, weights):
    c_prime = model.extract_capacitance_matrix_per_length(names4)
    c_node = c_prime * L_EL
    cap = CapacitanceMatrix(C_node=c_node, T=T)
    c_diff = cap.C_diff
    
    bc_x = {'+x': 0.5, '-x': -0.5, '+y': 0.0, '-y': 0.0}
    ux = model.solve_potential(bc_x)
    E_x = model.extract_electric_field(ux, pts)
    
    bc_y = {'+x': 0.0, '-x': 0.0, '+y': 0.5, '-y': -0.5}
    uy = model.solve_potential(bc_y)
    E_y = model.extract_electric_field(uy, pts)
    
    A = np.zeros((2, 2))
    A[0, 0] = np.sum(E_x[0] * weights)
    A[1, 0] = np.sum(E_x[1] * weights)
    A[0, 1] = np.sum(E_y[0] * weights)
    A[1, 1] = np.sum(E_y[1] * weights)
    
    rms = np.zeros((2, 2))
    rms[0, 0] = np.sqrt(np.sum((E_x[0] - A[0, 0])**2 * weights))
    rms[1, 0] = np.sqrt(np.sum((E_x[1] - A[1, 0])**2 * weights))
    rms[0, 1] = np.sqrt(np.sum((E_y[0] - A[0, 1])**2 * weights))
    rms[1, 1] = np.sqrt(np.sum((E_y[1] - A[1, 1])**2 * weights))
    
    vpi_x = EPI / abs(A[0, 0]) if A[0, 0] != 0 else float('inf')
    vpi_y = EPI / abs(A[1, 1]) if A[1, 1] != 0 else float('inf')
    
    return {
        'c_node': c_node, 'c_diff': c_diff,
        'A': A, 'vpi_x': vpi_x, 'vpi_y': vpi_y,
        'rms_xx_pct': abs(rms[0, 0] / A[0, 0]) * 100 if A[0, 0] != 0 else 0,
        'rms_yy_pct': abs(rms[1, 1] / A[1, 1]) * 100 if A[1, 1] != 0 else 0,
        'leak_xy': abs(A[0, 1] / A[0, 0]) * 100 if A[0, 0] != 0 else 0,
        'leak_yx': abs(A[1, 0] / A[1, 1]) * 100 if A[1, 1] != 0 else 0,
        'uy': uy
    }

def run_decision():
    setbacks = [500, 750]
    freqs = [6e6, 10e6, 50e6, 100e6]
    offsets = {
        'centered': (0, 0),
        '+100x': (100e-6, 0),
        '-100x': (-100e-6, 0),
        '+100y': (0, 100e-6),
        '-100y': (0, -100e-6),
        '+200x': (200e-6, 0),
        '-200x': (-200e-6, 0),
        '+200y': (0, 200e-6),
        '-200y': (0, -200e-6),
    }

    results = {}
    
    for sb in setbacks:
        print(f"Analyzing {sb} um setback...")
        strip_width = 3000 - 2 * sb
        mesh = build_hybrid_mesh(CS, sb * 1e-6, res=RES)
        model = ElectrostaticModel2D(mesh, lambda p: np.full(p.shape[1], EPS_R))
        
        # Center analysis
        pts0, w0 = make_eval_grid(0, 0)
        m0 = analyze_config(model, pts0, w0)
        
        rf_currents = {}
        for f in freqs:
            omega = 2 * np.pi * f
            # Quadrature drive currents: I = omega * |C_diff @ V_diff|
            # V_diff = [Vpi_x, -1j * Vpi_y]
            # Ix = j*omega*(Cxx*Vpi_x - 1j*Cxy*Vpi_y)
            # |Ix| = omega * sqrt((Cxx*Vpi_x)^2 + (Cxy*Vpi_y)^2)
            # Iy = j*omega*(Cyx*Vpi_x - 1j*Cyy*Vpi_y)
            # |Iy| = omega * sqrt((Cyx*Vpi_x)^2 + (Cyy*Vpi_y)^2)
            ix_quad = omega * np.sqrt((m0['c_diff'][0,0] * m0['vpi_x'])**2 + (m0['c_diff'][0,1] * m0['vpi_y'])**2)
            iy_quad = omega * np.sqrt((m0['c_diff'][1,0] * m0['vpi_x'])**2 + (m0['c_diff'][1,1] * m0['vpi_y'])**2)
            rf_currents[f] = {'ix': ix_quad, 'iy': iy_quad}
            
        offset_res = {}
        for name, (dx, dy) in offsets.items():
            pts, w = make_eval_grid(dx, dy)
            m = analyze_config(model, pts, w)
            offset_res[name] = m
            
        results[sb] = {
            'strip_width': strip_width,
            'center': m0,
            'rf': rf_currents,
            'offsets': offset_res,
            'model': model
        }
        
    return results

if __name__ == '__main__':
    res = run_decision()
    
    lines = []
    W = lines.append
    
    W("# 500 µm vs 750 µm Hybrid Geometry Decision Report")
    W("")
    W("Comparison of added-electrode setback options for the +y/-y strip electrodes.")
    W("")
    
    for sb in [500, 750]:
        r = res[sb]
        c = r['center']
        rf = r['rf']
        W(f"## Option: {sb} µm Setback")
        W(f"**1. Strip Width:** {r['strip_width']} µm")
        W("")
        W("**2. Cdiff Matrix (pF):**")
        W(f"- C_xx: {c['c_diff'][0,0]*1e12:.3f}")
        W(f"- C_yy: {c['c_diff'][1,1]*1e12:.3f}")
        W(f"- C_xy: {c['c_diff'][0,1]*1e12:.3f}")
        W("")
        W("**3. A Matrix (V/m/V):**")
        W(f"- A_xx: {c['A'][0,0]:.2f}")
        W(f"- A_yy: {c['A'][1,1]:.2f}")
        W(f"- A_xy (from Vy): {c['A'][0,1]:.2f}")
        W(f"- A_yx (from Vx): {c['A'][1,0]:.2f}")
        W("")
        W("**4. Half-Wave Voltages:**")
        W(f"- $V_{{\\pi,x}}$: {c['vpi_x']:.1f} V")
        W(f"- $V_{{\\pi,y}}$: {c['vpi_y']:.1f} V")
        W("")
        W(f"**5. Quadrature Amplitude Ratio:** $V_{{\\pi,y}} / V_{{\\pi,x}}$ = {c['vpi_y']/c['vpi_x']:.3f}")
        W("")
        W("**6. Quadrature Peak Currents (A):**")
        W("| Frequency (MHz) | $I_{peak,x}$ | $I_{peak,y}$ |")
        W("|-----------------|-------------|-------------|")
        for f in [6e6, 10e6, 50e6, 100e6]:
            W(f"| {int(f/1e6)} | {rf[f]['ix']:.3f} | {rf[f]['iy']:.3f} |")
        W("")
        W("**7 & 8. Nonuniformity and Leakage vs Beam Offset:**")
        W("| Offset | RMS$_x$ (%) | RMS$_y$ (%) | Leak$_{xy}$ (%) | Leak$_{yx}$ (%) |")
        W("|--------|------------|------------|--------------|--------------|")
        for off in ['centered', '+100x', '-100x', '+100y', '-100y', '+200x', '-200x', '+200y', '-200y']:
            m = r['offsets'][off]
            W(f"| {off} | {m['rms_xx_pct']:.2f} | {m['rms_yy_pct']:.2f} | {m['leak_xy']:.3f} | {m['leak_yx']:.3f} |")
        W("")
        W("---")
        W("")
        
    # Generate field maps for y-drive
    for sb in [500, 750]:
        model = res[sb]['model']
        uy = res[sb]['center']['uy']
        
        Nf = 100
        xf = np.linspace(-CS/2, CS/2, Nf)
        yf = np.linspace(-CS/2, CS/2, Nf)
        Xf, Yf = np.meshgrid(xf, yf)
        pts_f = np.vstack((Xf.ravel(), Yf.ravel()))
        
        E_ydrive = model.extract_electric_field(uy, pts_f)
        Ey_map = E_ydrive[1].reshape(Nf, Nf)
        
        fig, ax = plt.subplots(figsize=(6, 5))
        ext = [-HALF*1e3, HALF*1e3, -HALF*1e3, HALF*1e3]
        im = ax.imshow(Ey_map, extent=ext, origin='lower', cmap='RdBu_r', aspect='equal')
        ax.set_title(f'$E_y$ from y-drive (setback = {sb} µm)', fontsize=12)
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        plt.colorbar(im, ax=ax, label='V/m per $V_{diff}$')
        
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(W0*1e3*np.cos(theta), W0*1e3*np.sin(theta), 'k--', linewidth=1, alpha=0.5)
        
        sb_mm = sb / 1000
        # Full face left/right
        ax.plot([-1.5, -1.5], [-1.5, 1.5], 'k-', linewidth=3)
        ax.plot([1.5, 1.5], [-1.5, 1.5], 'k-', linewidth=3)
        # Strip top/bottom
        ax.plot([-1.5+sb_mm, 1.5-sb_mm], [1.5, 1.5], 'r-', linewidth=3)
        ax.plot([-1.5+sb_mm, 1.5-sb_mm], [-1.5, -1.5], 'r-', linewidth=3)
        
        fig.tight_layout()
        fig.savefig(out_dir / f"decision_fieldmap_y_s{sb}.png", dpi=150)
        
    plt.close('all')
    
    W("## Y-Drive Field Maps")
    W("![500 um](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/decision_fieldmap_y_s500.png)")
    W("![750 um](C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680/decision_fieldmap_y_s750.png)")
    W("")
    
    # Final recommendation logic
    W("## 9. Final Recommendation")
    W("")
    W("**Stricter first-build rules:**")
    W("- $V_{\\pi,y} < 750$ V")
    W("- RMS$_y < 7\\%$ for centered and $\\pm 100$ µm offsets")
    W("- RMS$_y < 10\\%$ for $\\pm 200$ µm offsets")
    W("- Cross-axis leakage $< 1\\%$")
    W("- Minimize max($I_x, I_y$) at 100 MHz")
    W("")
    
    W("| Rule | 500 µm | 750 µm |")
    W("|------|--------|--------|")
    
    def check_rules(r):
        vpi_pass = r['center']['vpi_y'] < 750
        rms7_pass = all(r['offsets'][k]['rms_yy_pct'] < 7 for k in ['centered', '+100x', '-100x', '+100y', '-100y'])
        rms10_pass = all(r['offsets'][k]['rms_yy_pct'] < 10 for k in ['+200x', '-200x', '+200y', '-200y'])
        leak_pass = all(max(r['offsets'][k]['leak_xy'], r['offsets'][k]['leak_yx']) < 1 for k in r['offsets'])
        return vpi_pass, rms7_pass, rms10_pass, leak_pass

    rules500 = check_rules(res[500])
    rules750 = check_rules(res[750])
    
    W(f"| $V_{{\\pi,y}} < 750$ V | {'PASS' if rules500[0] else 'FAIL'} | {'PASS' if rules750[0] else 'FAIL'} |")
    W(f"| RMS$_y < 7\\%$ (0/100) | {'PASS' if rules500[1] else 'FAIL'} | {'PASS' if rules750[1] else 'FAIL'} |")
    W(f"| RMS$_y < 10\\%$ (200) | {'PASS' if rules500[2] else 'FAIL'} | {'PASS' if rules750[2] else 'FAIL'} |")
    W(f"| Leakage $< 1\\%$ | {'PASS' if rules500[3] else 'FAIL'} | {'PASS' if rules750[3] else 'FAIL'} |")
    W("")
    
    if all(rules500) and all(rules750):
        ix5 = res[500]['rf'][100e6]['ix']
        iy5 = res[500]['rf'][100e6]['iy']
        max5 = max(ix5, iy5)
        ix7 = res[750]['rf'][100e6]['ix']
        iy7 = res[750]['rf'][100e6]['iy']
        max7 = max(ix7, iy7)
        if max5 < max7:
            W(f"**Selected:** 500 µm setback. Both passed, but 500 µm has lower max current ({max5:.2f} A vs {max7:.2f} A).")
        else:
            W(f"**Selected:** 750 µm setback. Both passed, but 750 µm has lower max current ({max7:.2f} A vs {max5:.2f} A).")
    elif all(rules500):
        W("**Selected:** **500 µm setback**. It passes all first-build safety margins, whereas 750 µm fails.")
    elif all(rules750):
        W("**Selected:** **750 µm setback**. It passes all first-build safety margins, whereas 500 µm fails.")
    else:
        W("**Warning:** Neither geometry passes all strict first-build rules.")
        
    with open(out_dir / "decision_report.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
