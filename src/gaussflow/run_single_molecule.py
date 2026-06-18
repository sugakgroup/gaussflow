import argparse
import json
import math
import sys
from pathlib import Path
from time import time
from rdkit import Chem

proj_root = (Path.cwd().parent).resolve()
sys.path.insert(0, str(proj_root))

from gaussflow.geometry import smiles_to_xyz_MMFF
from gaussflow.gaussian_singlerun import run_single
from gaussflow.gaussian_parser import output_to_property

def main():
    time0 = time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--id", type=int, required=True)

    idx = parser.parse_args().idx

    base_dir = Path(f"output/mol_{idx}")
    config = json.loads((base_dir / "config.json").read_text(encoding="utf-8"))
    properties = dict()
    properties["molecule"] = config["molecule"]
    out_file = base_dir / f"results_{idx}.out"

    for item in config["workflow"]:
        with open(out_file, "a") as f:
            f.write(f"Running {item['name']}...\n")
        if item["tool"] == "geometry":
            if item["method"] == "default_v0.1.0":
                success, results = smiles_to_xyz_MMFF(config["molecule"][item["structure_source"]])
                properties[item["name"]] = output_to_property(results, item["parse"]["target_properties"])

        elif item["tool"] == "g16":
            xyz = results
            for inner in item["structure_source"].split(":")[::-1]:
                xyz = results[inner]
            success, results = run_single(
                xyz=xyz,
                id=idx,
                wfname=item.get("name","test"),
                target_properties=item["parse"]["target_properties"],
                nprocs=config["compute"]["cpus_per_task"],
                mem=config["compute"]["mem"],
                route_line=item["route_line"],
                charge=item.get("charge", 0),
                multiplicity=item.get("multiplicity", 1),
                remove_results=False
            )
            properties[item["name"]] = output_to_property(results=results, target_properties=item["parse"]["target_properties"])
        else:
            raise ValueError(f"Unknown tool: {item['tool']}")

        with open(out_file, "a") as f:
            f.write(f"Finished {item['name']}.\n")
            f.write(f"Success: {success}\n")
            f.write(f"Elapsed time: {time() - time0:.2f} seconds\n")
            f.write(f"Results for {item['name']}: {results}\n")
            f.write(f"Properties for {item['name']}: {properties[item['name']]}\n")
            if not success and item.get("stop_on_fail", False):
                f.write(f"Workflow stopped due to failure in {item['name']}.")
                break
            else:
                f.write(f"----------------------------------------------------------------\n")
        
    # save results
    with open(base_dir / f"properties.json", "w") as f:
        json.dump(properties, f, indent=4)

if __name__ == "__main__":
    main()