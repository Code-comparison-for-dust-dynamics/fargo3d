#!/bin/bash
#SBATCH --job-name=mB-dust
#SBATCH --account=def-evelee
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out
# Submit from the repository root on trig-login01. See cluster/README.md.
source "${SLURM_SUBMIT_DIR:?submit from the repository root}/cluster/dust_body.sh" \
  "${SLURM_SUBMIT_DIR}/cluster/profiles/trillium.sh"
