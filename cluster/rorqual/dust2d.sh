#!/bin/bash
#SBATCH --job-name=mB2D-dust
#SBATCH --account=def-evelee
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus-per-node=h100:1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x-%j.out
# 2D (R, phi) Model B on one GPU. Submit from the repository root. See cluster/README.md.
export DIM=2D
source "${SLURM_SUBMIT_DIR:?submit from the repository root}/cluster/dust_body.sh" \
  "${SLURM_SUBMIT_DIR}/cluster/profiles/rorqual.sh"
