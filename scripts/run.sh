#!/bin/bash
#SBATCH --job-name=suga
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=/work/log/slurm-%j.out
#SBATCH --error=/work/log/slurm-%j.err

source /home/suga/miniconda3/etc/profile.d/conda.sh
conda activate CARBOT_base

pip install -e .
gaussflow submit --config config/workflow_T0_S1_T1.json