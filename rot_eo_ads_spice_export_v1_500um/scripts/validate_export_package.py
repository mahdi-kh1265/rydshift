import os
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
