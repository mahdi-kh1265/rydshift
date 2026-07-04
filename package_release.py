import os
import sys
import shutil
import hashlib
from datetime import datetime
from pathlib import Path

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def main():
    base_dir = Path("c:/Users/khams008/Documents/rydshift/rot_eo_jupyter_handoff")
    release_dir = base_dir / "rot_eo_ads_spice_export_v1_500um"
    
    if release_dir.exists():
        shutil.rmtree(release_dir)
        
    # 1. Create directory structure
    dirs = [
        "geometry_source",
        "ADS",
        "SPICE_LTspice_ngspice",
        "validation",
        "scripts",
        "checksums"
    ]
    for d in dirs:
        (release_dir / d).mkdir(parents=True, exist_ok=True)
        
    # 2. Geometry Source
    shutil.copy2(base_dir / "rf_ready_geometry_spec.json", release_dir / "geometry_source" / "rf_ready_geometry_spec.json")
    # Write a simple markdown summary and port map
    with open(release_dir / "geometry_source" / "geometry_summary.md", "w", encoding="utf-8") as f:
        f.write("# Geometry Summary: Locked 500 µm Hybrid Setback\n\n")
        f.write("- **Material**: z-cut LiNbO3\n")
        f.write("- **Dimensions**: 3 mm × 3 mm × 30 mm\n")
        f.write("- **x-electrodes**: Full-face left and right (width = 3 mm)\n")
        f.write("- **y-electrodes**: Centered top and bottom strips (width = 2 mm, setback = 500 µm)\n")
        f.write("- **Optical Beam**: 500 µm Gaussian waist\n")
        
    with open(release_dir / "geometry_source" / "port_map.md", "w", encoding="utf-8") as f:
        f.write("# Port Definitions\n\n")
        f.write("## 4-Port Single-Ended\n")
        f.write("1. `+x_full` (left face)\n")
        f.write("2. `-x_full` (right face)\n")
        f.write("3. `+y_strip` (top face)\n")
        f.write("4. `-y_strip` (bottom face)\n\n")
        f.write("## 2-Port Differential\n")
        f.write("1. `x` differential pair\n")
        f.write("2. `y` differential pair\n")

    # 3. ADS
    ads_src = base_dir / "ads_export"
    ads_dest = release_dir / "ADS"
    # Copy relevant Touchstone and wrappers
    for f in ads_src.glob("*.s4p"): shutil.copy2(f, ads_dest)
    for f in ads_src.glob("*.s2p"): shutil.copy2(f, ads_dest)
    for f in ads_src.glob("*.s1p"): shutil.copy2(f, ads_dest)
    shutil.copy2(ads_src / "fixture_parasitics_wrapper.net", ads_dest)
    
    with open(ads_dest / "ADS_FIRST_USE_CHECKLIST.md", "w", encoding="utf-8") as f:
        f.write("# ADS First Use Checklist\n\n")
        f.write("Primary model: `eo_load_4port_single_ended_R50.s4p`\n")
        f.write("Port 1 = +x_full\n")
        f.write("Port 2 = -x_full\n")
        f.write("Port 3 = +y_strip\n")
        f.write("Port 4 = -y_strip\n\n")
        f.write("- [ ] Place an SnP component pointing to `eo_load_4port_single_ended_R50.s4p`.\n")
        f.write("- [ ] Ground or terminate ports with 50 ohms to verify low-frequency capacitance.\n")
        f.write("- [ ] **Status**: Pending actual simulator import.\n")

    # 4. SPICE
    spice_dest = release_dir / "SPICE_LTspice_ngspice"
    for f in ads_src.glob("*.sp"): shutil.copy2(f, spice_dest)
    for f in ads_src.glob("*.sub"): shutil.copy2(f, spice_dest)
    for f in ads_src.glob("*.cir"): shutil.copy2(f, spice_dest)
    
    with open(spice_dest / "SPICE_FIRST_USE_CHECKLIST.md", "w", encoding="utf-8") as f:
        f.write("# SPICE / LTspice / ngspice First Use Checklist\n\n")
        f.write("- **Pure capacitor model** (`eo_load_4node.sp`) is the portable baseline.\n")
        f.write("- **Fixed resistor lossy models** (`...at_100MHz.sp`) are exact **only** at 100 MHz.\n")
        f.write("- **RC-ladder lossy models** (`...ladder_0p5MHz_200MHz.sub`) are finite-band approximations over 0.5–200 MHz, great for transient sweeps.\n")
        f.write("- To instantiate: Use `.INCLUDE \"eo_load_4node.sp\"` and instance it as `X1 p_x m_x p_y m_y EO_LOAD_4NODE`.\n")
        f.write("- **LTspice Status**: Pending actual simulator import.\n")

    # 5. Scripts
    scripts_dest = release_dir / "scripts"
    shutil.copy2(ads_src / "generate_ads_exports.py", scripts_dest)
    shutil.copy2(ads_src / "fit_vna_measurements.py", scripts_dest)
    
    # Create validate_export_package.py using relative paths
    validation_script_content = """import os
import json
import skrf as rf
from pathlib import Path

def main():
    print("Running portable validation inside the release folder...")
    base_dir = Path(__file__).parent.parent
    ads_dir = base_dir / "ADS"
    
    # Just a simple sanity check script that can run independently
    s4p_file = ads_dir / "eo_load_4port_single_ended_R50.s4p"
    if s4p_file.exists():
        ntwk = rf.Network(str(s4p_file))
        print(f"PASS: Read {s4p_file.name} successfully (Ports: {ntwk.number_of_ports}, Freq: {ntwk.f[0]/1e6} - {ntwk.f[-1]/1e6} MHz)")
    else:
        print("FAIL: Core s4p file missing!")
        
if __name__ == '__main__':
    main()
"""
    with open(scripts_dest / "validate_export_package.py", "w", encoding="utf-8") as f:
        f.write(validation_script_content)

    # 6. Validation Reports
    val_dest = release_dir / "validation"
    shutil.copy2(ads_src / "validation_plots.png", val_dest)
    
    # We will copy the validation reports from the artifacts folder
    artifact_dir = Path("C:/Users/khams008/.gemini/antigravity-ide/brain/16a39c4f-cd79-46a8-9129-286bf7ecb680")
    for r in ["ads_export_strict_validation.md", "extended_audit_report.md"]:
        if (artifact_dir / r).exists():
            shutil.copy2(artifact_dir / r, val_dest)
            
    # Also grab the compatibility report we just made
    if (base_dir / "ads_spice_format_compatibility_audit.md").exists():
        shutil.copy2(base_dir / "ads_spice_format_compatibility_audit.md", val_dest)

    # 7. Write README_START_HERE.md
    with open(release_dir / "README_START_HERE.md", "w", encoding="utf-8") as f:
        f.write("# rot_eo_ads_spice_export_v1_500um\n\n")
        f.write("## What is this package?\n")
        f.write("This is a curated, versioned export of the RF-ready linear load models for the rotating electro-optic (EO) frequency shifter. It contains machine-readable Touchstone (.sNp) and SPICE (.sp/.sub) files.\n\n")
        f.write("## What geometry does it represent?\n")
        f.write("It represents the locked 500 µm hybrid setback geometry (z-cut LiNbO3, 3x3x30 mm).\n\n")
        f.write("## Which file should I use first in ADS?\n")
        f.write("Use `ADS/eo_load_4port_single_ended_R50.s4p` using the SnP component. See `ADS_FIRST_USE_CHECKLIST.md`.\n\n")
        f.write("## Which file should I use first in LTspice/ngspice?\n")
        f.write("Use `SPICE_LTspice_ngspice/eo_load_4node.sp`. See `SPICE_FIRST_USE_CHECKLIST.md`.\n\n")
        f.write("## Port Definitions\n")
        f.write("See `geometry_source/port_map.md` for definitions.\n\n")
        f.write("## Differences between models\n")
        f.write("- **Lossless**: Pure authoritative capacitance matrix.\n")
        f.write("- **Lossy (RC-Ladder)**: Passive broadband approximation of a constant loss tangent.\n")
        f.write("- **Tolerance**: Scaled bounds (±20%) for extreme operating corner checking.\n")
        f.write("- **Fixture Wrapper**: A template for appending physical parasitics (wirebonds/traces).\n\n")
        f.write("## Validation Status\n")
        f.write("- **Math & Formatting**: Validated perfectly against Touchstone/SPICE specifications.\n")
        f.write("- **ADS/LTspice Smoke Testing**: **PENDING** actual simulator import by the RF team.\n")

    # 8. MANIFEST.md
    with open(release_dir / "MANIFEST.md", "w", encoding="utf-8") as f:
        f.write("# Manifest\n\n")
        for root, dirs, files in os.walk(release_dir):
            for file in files:
                filepath = Path(root) / file
                rel_path = filepath.relative_to(release_dir)
                if rel_path.name in ["MANIFEST.md", "README_START_HERE.md", "VERSION.txt"]:
                    f.write(f"- `{rel_path}`: Package documentation\n")
                elif "s4p" in rel_path.name:
                    f.write(f"- `{rel_path}`: Touchstone (ADS/RF), Authoritative or derived\n")
                elif "sp" in rel_path.name or "sub" in rel_path.name:
                    f.write(f"- `{rel_path}`: SPICE netlist, Authoritative or derived lossy\n")
                else:
                    f.write(f"- `{rel_path}`: Supporting asset/script\n")
                    
    # 9. VERSION.txt
    generator_hash = compute_sha256(ads_src / "generate_ads_exports.py")
    with open(release_dir / "VERSION.txt", "w", encoding="utf-8") as f:
        f.write("Package: rot_eo_ads_spice_export\n")
        f.write("Version: v1_500um\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n")
        f.write("Geometry: Locked hybrid 500 µm setback\n")
        f.write("Source: rf_ready_geometry_spec.json\n")
        f.write(f"Generator Hash: {generator_hash}\n")
        
    # 10. Checksums
    checksum_file = release_dir / "checksums" / "SHA256SUMS.txt"
    with open(checksum_file, "w", encoding="utf-8") as f:
        for root, dirs, files in os.walk(release_dir):
            for file in sorted(files):
                if file == "SHA256SUMS.txt": continue
                filepath = Path(root) / file
                rel_path = filepath.relative_to(release_dir)
                # Only checksum specific machine-readable files and scripts
                if filepath.suffix in [".s1p", ".s2p", ".s4p", ".sp", ".sub", ".net", ".cir", ".json", ".py", ".md", ".txt"]:
                    f.write(f"{compute_sha256(filepath)}  {rel_path.as_posix()}\n")
                    
    print(f"Created release folder: {release_dir.name}")

if __name__ == '__main__':
    main()
