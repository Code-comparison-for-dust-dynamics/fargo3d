#!/bin/bash
#SBATCH --job-name=mB2D-gas
#SBATCH --account=def-evelee
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=h100:1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=slurm-%x-%j.out
# 2D (R, phi) Model B on one GPU. Submit from the repository root. See cluster/README.md.
export DIM=2D
source "${SLURM_SUBMIT_DIR:?submit from the repository root}/cluster/gas_body.sh" \
  "${SLURM_SUBMIT_DIR}/cluster/profiles/nibi.sh"
