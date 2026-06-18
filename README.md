# gaussflow

A lightweight Gaussian/Slurm workflow manager developed for SugaGroup and reusable high-throughput quantum chemistry calculations.

`gaussflow` is a small, config-driven Python package for preparing and submitting Gaussian workflows on a Slurm cluster. It is currently designed around the SugaGroup computing environment, but the workflow structure is intentionally simple so that it can be adapted to other Gaussian/Slurm setups.

> Status: early pre-release. APIs, configuration keys, directory layout, and scheduler assumptions may change.

## What gaussflow does

`gaussflow` automates a minimal calculation loop:

1. Read a molecule table from CSV.
2. Create one working directory per molecule.
3. Generate a per-molecule `config.json`.
4. Submit each molecule workflow to Slurm.
5. Run workflow steps such as geometry generation and Gaussian calculations.
6. Parse selected Gaussian outputs.
7. Submit a dependent post-processing job.
8. Collect selected properties into `summary.csv`.

## Features

* CSV-based molecule input
* Config-driven workflow definition
* RDKit/MMFF-based initial 3D geometry generation from SMILES
* Gaussian 16 job execution
* Slurm job submission
* Slurm dependency job for after-run parsing
* Per-molecule `properties.json` output
* Final `summary.csv` generation

## Installation

```bash
git clone https://github.com/sugakgroup/gaussflow.git
cd gaussflow
pip install -e .
```

Python 3.10 or later is required.

## External requirements

* Gaussian 16
* Slurm
* RDKit
* A Python environment available on compute nodes

The current scheduler script assumes the SugaGroup cluster environment, including the conda activation command and Slurm defaults. For use on another cluster, edit `src/gaussflow/scheduler.py`.

## Basic usage

```bash
gaussflow submit --config config.json
```

## Example molecule table

```csv
id,smiles
0,C1=CC=CC=C1
1,C1=CC=C2C=CC=CC2=C1
```

## Example config

```json
{
  "project": {
    "name": "workflow_T0_S1_T1"
  },
  "input": {
    "molecule_table": "input/molecules.csv",
    "id_column": "id"
  },
  "output": {
    "root": "output"
  },
  "compute": {
    "cpus_per_task": 8,
    "mem": "8GB"
  },
  "workflow": [
    {
      "name": "geom",
      "tool": "geometry",
      "method": "default_v0.1.0",
      "structure_source": "smiles"
      "stop_on_fail": true
    },
    {
      "name": "opt",
      "tool": "g16",
      "structure_source": "xyz:geom",
      "route_line": "#p b3lyp/6-31g(d) opt",
      "charge": 0,
      "multiplicity": 1,
      "parse": {
        "target_properties": ["xyz", "SCF_energy"]
      },
      "stop_on_fail": true
    }
  ],
  "after_jobs": {
    "summary_columns": [
      "id:molecule",
      "smiles:molecule",
      "SCF_energy:opt"
    ]
  }
}
```

## Output layout

```text
output/
  mol_0/
    opt/
      mol_0.gjf
      mol_0.chk
      mol_0.log
    0.log
    config.json
    properties.json
    submit.sh
    results_0.out
  mol_1/
    opt/
      mol_1.gjf
      mol_1.chk
      mol_1.log
    1.log
    config.json
    properties.json
    submit.sh
    results_1.out
  after_jobs.sh
  summary.csv
```

## Development status

This is an early SugaGroup-oriented pre-release. The code is useful as a lightweight starting point for Gaussian workflow automation, but the workflow schema, CLI, parser interface, and scheduler abstraction are expected to evolve.

## License

This project is licensed under the MIT License.
