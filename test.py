from gaussflow.gaussian_parser import analyze_g16_output
from pathlib import Path

if __name__ == "__main__":
    print(analyze_g16_output(id=0, base_dir=Path("output/mol_0"), wfname="s1t1sp", target_properties=["xyz","SCF_energy","s1-t1_gap","excitation_energy_s1","excitation_energy_t1"]))