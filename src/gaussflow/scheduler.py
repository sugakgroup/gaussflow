# src/gaussflow/scheduler.py

from __future__ import annotations

from pathlib import Path
import subprocess


def write_slurm_script(output_root, id, config):
    script_path = output_root / f"mol_{id}" / "submit.sh"
    log_path = output_root / f"mol_{id}" / f"{id}.log"

    partition = config.get("scheduler", {}).get("partition", "cpu")
    nodelist = config.get("compute", {}).get("nodelist", "suga[01-12]")
    nodes = config.get("compute", {}).get("nodes", 1)
    ntasks = config.get("compute", {}).get("ntasks", 1)
    cpus_per_task = config.get("compute", {}).get("cpus_per_task", 1)
    mem = config.get("compute", {}).get("mem", "1GB")
    
    script = f"""#!/bin/bash
        #SBATCH --job-name=mol_{id}
        #SBATCH --partition={partition}
        #SBATCH --nodelist={nodelist}
        #SBATCH --nodes={nodes}
        #SBATCH --ntasks={ntasks}
        #SBATCH --cpus-per-task={cpus_per_task}
        #SBATCH --mem={mem}
        #SBATCH --output=/work/log/slurm-%j.out
        #SBATCH --error=/work/log/slurm-%j.err
        source /home/suga/miniconda3/etc/profile.d/conda.sh
        conda activate CARBOT_base
        
        python gaussflow.run_single_molecule.py --id {id} > {log_path} 2>&1
        """

    script_path.write_text(script, encoding="utf-8")
    script_path.chmod(0o755)

    return script_path


def submit_slurm_job(script_path):
    result = subprocess.run(
        ["sbatch", str(script_path)],
        check=True,
        text=True,
        capture_output=True,
    )

    # Usually: "Submitted batch job 123456"
    return result.stdout.strip().split()[-1]


def write_after_jobs_script(output_root, config_path):
    script_path = output_root / "after_jobs.sh"
    script = f"""#!/bin/bash
        #SBATCH --job-name=after_job
        #SBATCH --partition=cpu
        #SBATCH --nodelist=suga[01-12]
        #SBATCH --nodes=1
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=10GB
        #SBATCH --output=/work/log/slurm-%j.out
        #SBATCH --error=/work/log/slurm-%j.err
        source /home/suga/miniconda3/etc/profile.d/conda.sh
        conda activate CARBOT_base
        
        python gaussflow.run_parse_after_jobs.py --config {config_path}
        """

    script_path.write_text(script, encoding="utf-8")
    script_path.chmod(0o755)

    return script_path

def submit_after_jobs(output_root, dependency_job_ids, config_path, mode="afterok"):
    """
    Submit a Slurm job that starts after all dependency jobs finish.

    mode:
        afterany: run after all jobs finish, regardless of success/failure
        afterok:  run only if all jobs finish successfully
    """
    if not dependency_job_ids:
        raise ValueError("dependency_job_ids is empty.")

    dependency = f"{mode}:" + ":".join(dependency_job_ids)

    script_path = write_after_jobs_script(output_root, config_path)

    result = subprocess.run(
        [
            "sbatch",
            f"--dependency={dependency}",
            str(script_path),
        ],
        check=True,
        text=True,
        capture_output=True,
    )

    job_id = result.stdout.strip().split()[-1]
    return job_id