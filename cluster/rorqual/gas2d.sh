#!/bin/bash
#SBATCH --job-name=mB2D-gas
#SBATCH --account=def-evelee
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus=h100_1g.10gb:1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x-%j.out
# 2D (R, phi) Model B on one GPU. Submit from the repository root. See cluster/README.md.
export DIM=2D
source "${SLURM_SUBMIT_DIR:?submit from the repository root}/cluster/gas_body.sh" \
  "${SLURM_SUBMIT_DIR}/cluster/profiles/rorqual.sh"
