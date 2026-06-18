# src/gaussflow/workflow.py
import os
from pathlib import Path

from gaussflow.config import load_config, dump_config
from gaussflow.molecule import load_molecules
from gaussflow.gaussian_input import write_gaussian_input
from gaussflow.scheduler import submit_after_jobs, write_slurm_script, submit_slurm_job


def submit_workflow(config_path):
    config = load_config(config_path)
    molecules = load_molecules(config)

    output_root = Path(config["output"]["root"])

    job_ids = []

    # submit each molecule's workflow
    for mol in molecules:
        mol_config = config.copy()
        mol_config["molecule"] = mol
        os.makedirs(output_root / f"mol_{mol['id']}", exist_ok=True)
        dump_config(config=mol_config, config_path=output_root / f"mol_{mol['id']}" / "config.json")

        submit_script = write_slurm_script(
            output_root=output_root,
            id=mol["id"],
            config=mol_config
        )
        job_id = submit_slurm_job(submit_script)
        job_ids.append(job_id)
        print(f"Submitted mol_{mol['id']}: job_id={job_id}")
    
    # summarize data
    parse_job_id = submit_after_jobs(
        output_root=output_root,
        dependency_job_ids=job_ids,
        config_path=config_path,
        mode="afterany",
    )


            