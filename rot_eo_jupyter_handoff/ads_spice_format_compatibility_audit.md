# ADS/SPICE Format Compatibility Audit

## 1. Touchstone Syntax Audit
### eo_load_xdiff_R50.s1p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (1) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 0.000e+00
- PASS: Max singular value: 1.000000 (Passive: True)

### eo_load_ydiff_R50.s1p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (1) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 0.000e+00
- PASS: Max singular value: 1.000000 (Passive: True)

### eo_load_2port_diff_R100.s2p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 100.0 `
- PASS: Port count (2) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 2.199e-16
- PASS: Max singular value: 1.000000 (Passive: True)

### eo_load_2port_diff_R50.s2p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (2) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 2.199e-16
- PASS: Max singular value: 1.000000 (Passive: True)

### eo_load_2port_diff_R50_tand_0.0001.s2p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (2) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 2.199e-16
- PASS: Max singular value: 0.999999 (Passive: True)

### eo_load_2port_diff_R50_tand_0.0005.s2p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (2) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 2.198e-16
- PASS: Max singular value: 0.999997 (Passive: True)

### eo_load_2port_diff_R50_tand_0.001.s2p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (2) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 2.197e-16
- PASS: Max singular value: 0.999995 (Passive: True)

### eo_load_4port_single_ended_R50.s4p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (4) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 5.389e-16
- PASS: Max singular value: 1.000000 (Passive: True)
- PASS: Max C_node reconstruction error: 6.291e-14 %

### eo_load_4port_single_ended_R50_tand_0.0001.s4p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (4) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 6.087e-16
- PASS: Max singular value: 1.000000 (Passive: True)
- PASS: Max C_node reconstruction error: 6.291e-14 %

### eo_load_4port_single_ended_R50_tand_0.0005.s4p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (4) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 5.551e-16
- PASS: Max singular value: 1.000000 (Passive: True)
- PASS: Max C_node reconstruction error: 6.291e-14 %

### eo_load_4port_single_ended_R50_tand_0.001.s4p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (4) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 5.979e-16
- PASS: Max singular value: 1.000000 (Passive: True)
- PASS: Max C_node reconstruction error: 6.291e-14 %

### eo_load_4port_single_ended_R50_tol_0.8.s4p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (4) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 5.185e-16
- PASS: Max singular value: 1.000000 (Passive: True)

### eo_load_4port_single_ended_R50_tol_1.2.s4p
- PASS: Pure ASCII.
- PASS: Header option line valid: `# Hz S RI R 50.0 `
- PASS: Port count (4) matches extension.
- PASS: Frequencies strictly increasing.
- PASS: Reciprocity error: 5.525e-16
- PASS: Max singular value: 1.000000 (Passive: True)

## 2. SPICE Syntax Audit
### eo_load_4node.sp
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.
- PASS: SPICE netlist reconstructed C_node max error: 5.990e-19 F

### eo_load_4node_lossy_tand_0.0001.sp
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### eo_load_4node_lossy_tand_0.0001_at_100MHz.sp
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### eo_load_4node_lossy_tand_0.0005.sp
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### eo_load_4node_lossy_tand_0.0005_at_100MHz.sp
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### eo_load_4node_lossy_tand_0.001.sp
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### eo_load_4node_lossy_tand_0.001_at_100MHz.sp
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### eo_load_4node_lossy_tand_0.0001_ladder_0p5MHz_200MHz.sub
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### eo_load_4node_lossy_tand_0.0005_ladder_0p5MHz_200MHz.sub
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### eo_load_4node_lossy_tand_0.001_ladder_0p5MHz_200MHz.sub
- PASS: Pure ASCII.
- PASS: .SUBCKT EO_LOAD_4NODE and .ENDS verified.

### fixture_parasitics_wrapper.net
- PASS: Pure ASCII.

### test_ltspice_core_cap.cir
- PASS: Pure ASCII.

### test_ngspice_core_cap.cir
- PASS: Pure ASCII.

### test_ngspice_ladder_loss.cir
- PASS: Pure ASCII.

## 3. Final Compatibility Status Table
| File | Math Validated | Syntax Validated | skrf/ngspice Readback | ADS Import Smoke Test | LTspice Smoke Test |
|---|---|---|---|---|---|
| eo_load_xdiff_R50.s1p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_ydiff_R50.s1p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_2port_diff_R100.s2p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_2port_diff_R50.s2p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_2port_diff_R50_tand_0.0001.s2p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_2port_diff_R50_tand_0.0005.s2p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_2port_diff_R50_tand_0.001.s2p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_4port_single_ended_R50.s4p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_4port_single_ended_R50_tand_0.0001.s4p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_4port_single_ended_R50_tand_0.0005.s4p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_4port_single_ended_R50_tand_0.001.s4p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_4port_single_ended_R50_tol_0.8.s4p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_4port_single_ended_R50_tol_1.2.s4p | ✅ | ✅ | ✅ | Pending | N/A |
| eo_load_4node.sp | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.0001.sp | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.0001_at_100MHz.sp | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.0005.sp | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.0005_at_100MHz.sp | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.001.sp | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.001_at_100MHz.sp | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.0001_ladder_0p5MHz_200MHz.sub | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.0005_ladder_0p5MHz_200MHz.sub | ✅ | ✅ | ✅ | Pending | Pending |
| eo_load_4node_lossy_tand_0.001_ladder_0p5MHz_200MHz.sub | ✅ | ✅ | ✅ | Pending | Pending |
| fixture_parasitics_wrapper.net | ✅ | ✅ | ✅ | Pending | N/A |
| test_ltspice_core_cap.cir | ✅ | ✅ | ✅ | Pending | Pending |
| test_ngspice_core_cap.cir | ✅ | ✅ | ✅ | Pending | Pending |
| test_ngspice_ladder_loss.cir | ✅ | ✅ | ✅ | Pending | Pending |