#!/bin/bash
#SBATCH --job-name=mB-gas
#SBATCH --account=def-evelee
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --gpus-per-node=h100:4
#SBATCH --mem=160G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%j.out
# Submit from the repository root. See cluster/README.md.
source "${SLURM_SUBMIT_DIR:?submit from the repository root}/cluster/gas_body.sh" \
  "${SLURM_SUBMIT_DIR}/cluster/profiles/fir.sh"
