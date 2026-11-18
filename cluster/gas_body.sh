#!/bin/bash
# Gas stage of the Model B runs: 0 to 100 orbits (outputs 0 to 20), shared by
# every dust run. Usage (from a Slurm wrapper): source gas_body.sh PROFILE
# Output: $FARGO_OUTPUT_ROOT/modelB_gas (default $SCRATCH/fargo3d_rerun/modelB_gas).

profile=${1:?usage: source gas_body.sh PROFILE}
source "$profile"
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

setup=modelB_gas
out="$output_root/$setup"

if [[ -d "$out" ]]; then
  if [[ -f "$out/gasdens0.dat" && -n "$(last_complete_output "$out" 20 "$out/gasdens0.dat")" ]]; then
    echo "Gas stage already complete in $out; nothing to do."
    exit 0
  fi
  echo "ERROR: $out exists but holds no complete output 20; remove it to rerun." >&2
  exit 4
fi
mkdir -p "$out"

load_modules
bdir=$(build_setup "$setup" "$out")
write_provenance "$out" "$bdir"

cd "$bdir"
echo "=== gas stage, job $SLURM_JOB_ID, $(date --iso-8601=seconds)" | tee -a "$out/launch.log"
mpirun -np "$SLURM_NTASKS" ./fargo3d -o "OutputDir,@$out" "setups/$setup/$setup.par" \
  2>&1 | tee -a "$out/gas_run.log"

if [[ -z "$(last_complete_output "$out" 20 "$out/gasdens0.dat")" ]]; then
  echo "ERROR: the gas stage ended without a complete output 20 in $out." >&2
  exit 5
fi
echo "=== gas stage complete, $(date --iso-8601=seconds)" | tee -a "$out/launch.log"
