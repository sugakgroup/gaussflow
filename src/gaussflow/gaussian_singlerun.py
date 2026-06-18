import json
from pathlib import Path
import shutil
import subprocess
import os

from gaussflow.gaussian_input import make_g16_input
from gaussflow.gaussian_parser import analyze_g16_output

def run_single(xyz, id=None, wfname="test", base_dir=None, target_properties=["xyz","energy","pop"], nprocs=1, mem="4GB", route_line="#", charge=0, multiplicity=1, remove_results=False):

    success = False
    results = {prop: None for prop in target_properties}

    if base_dir is None:
        base_dir = Path(f"output/mol_{id}")

    if id is None:
        id = os.getpid()
    
    if len(xyz) == 0:
        return success, results
    os.makedirs(base_dir / wfname, exist_ok=True)

    make_g16_input(xyz,id,base_dir=base_dir,wfname=wfname,nprocs=nprocs,mem=mem,charge=charge,multiplicity=multiplicity,route_line=route_line)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(nprocs)
    proc = subprocess.Popen(["g16", f'mol_{id}.gjf', f'mol_{id}.log'], cwd=base_dir / wfname, stderr=subprocess.DEVNULL, env=env)
    proc.wait()

    success, results = analyze_g16_output(id,base_dir,wfname,target_properties)
    if not success:
        return success, results
    
    if remove_results:
        os.remove(base_dir / wfname / f'mol_{id}.gjf')
        os.remove(base_dir / wfname / f'mol_{id}.log')
        shutil.rmtree(base_dir / wfname)

    return success, results