import os

def make_g16_input(xyz,id,base_dir,wfname="test",nprocs=1,mem="4GB",charge=0,multiplicity=1,route_line="#"):
    os.makedirs(base_dir / wfname, exist_ok=True)
    with open(base_dir / wfname / f"mol_{id}.gjf", "w") as f:
        f.write(f"%nprocshared={nprocs}\n")
        f.write(f"%mem={mem}\n")
        f.write(f"%chk=mol_{id}.chk\n")
        f.write(f"{route_line}\n\n")
        f.write(f"{id}\n\n")
        f.write(f"{charge} {multiplicity}\n")
        for atom in xyz:
            f.write(f" {atom[0]} {atom[1]} {atom[2]} {atom[3]}\n")
        f.write(f"\n\n")