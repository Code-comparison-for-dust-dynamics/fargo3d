#!/bin/bash
# Dust stage of one Model B run, 100 to 1500 orbits (outputs 20 to 300), in
# segments: the first segment copies output 20 of the gas stage and releases
# the dust (-p); every later segment restarts from the last complete output.
# A segment that finds output 300 exits at once, so segments can be chained
# with --dependency=afterany. Usage (from a Slurm wrapper, with CASE and ST
# exported): source dust_body.sh PROFILE
#   CASE = contdiff | nodiff,  ST = 0.1 | 0.0316 | 0.01
#   GAS_DIR (optional) = gas-stage outputs, default $FARGO_OUTPUT_ROOT/modelB_gas
# Output: $FARGO_OUTPUT_ROOT/<CASE>_St<ST> (default root $SCRATCH/fargo3d_rerun); with DIM=2D, the
# setups modelB2D_* and outputs 2D_<CASE>_St<ST> (gas stage 2D_modelB_gas).

profile=${1:?usage: source dust_body.sh PROFILE}
source "$profile"
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

case "${CASE:?export CASE=contdiff or nodiff}" in
  contdiff|nodiff) ;;
  *) echo "ERROR: CASE must be contdiff or nodiff, got $CASE" >&2; exit 2 ;;
esac
case "${ST:?export ST=0.1, 0.0316 or 0.01}" in
  0.1) invstokes=10 ;;
  0.0316) invstokes=31.62 ;;
  0.01) invstokes=100 ;;
  *) echo "ERROR: ST must be 0.1, 0.0316 or 0.01, got $ST" >&2; exit 2 ;;
esac

setup=${SETUP_PREFIX}_$CASE
gas_dir=${GAS_DIR:-"$output_root/${OUT_PREFIX}modelB_gas"}
out="$output_root/${OUT_PREFIX}${CASE}_St$ST"
final=300

if [[ -f "$out/summary$final.dat" ]]; then
  echo "$out already reached output $final; nothing to do."
  exit 0
fi

if [[ -d "$out" && -f "$out/gasdens20.dat" ]]; then
  restart=$(last_complete_output "$out" 20 "$out/gasdens20.dat")
else
  restart=""
fi

if [[ -n "$restart" && -f "$out/provenance.txt" ]]; then
  first_commit=$(sed -n 's/^git_commit=//p' "$out/provenance.txt" | head -1)
  now_commit=$(git -C "$repo_root" rev-parse HEAD)
  if [[ "$first_commit" != "$now_commit" && -z "${ALLOW_COMMIT_CHANGE:-}" ]]; then
    echo "ERROR: $out was started at commit $first_commit but the repository is at $now_commit." >&2
    echo "       Check out $first_commit, or set ALLOW_COMMIT_CHANGE=1 if the change is intended." >&2
    exit 6
  fi
fi

if [[ -z "$restart" ]]; then
  # First segment: start from output 20 of the gas stage.
  if [[ -z "$(last_complete_output "$gas_dir" 20 "$gas_dir/gasdens20.dat" 2>/dev/null)" ]]; then
    echo "ERROR: no complete gas-stage output 20 in $gas_dir." >&2
    exit 4
  fi
  if [[ -d "$out" ]] && compgen -G "$out/*.dat" > /dev/null; then
    echo "ERROR: $out holds data but no complete output; remove it to start over." >&2
    exit 4
  fi
  mkdir -p "$out"
  cp "$gas_dir"/*20.dat "$gas_dir"/planet0.dat "$gas_dir"/domain_*.dat "$out"/
  rm -f "$out"/*_2d.dat
fi

load_modules
bdir=$(build_setup "$setup" "$out")
write_provenance "$out" "$bdir"

cd "$bdir"
opts="Invstokes1,$invstokes,OutputDir,@$out"
if [[ -z "$restart" ]]; then
  echo "=== $CASE St=$ST: dust release from output 20, job $SLURM_JOB_ID, $(date --iso-8601=seconds)" | tee -a "$out/launch.log"
  mpirun -np "$SLURM_NTASKS" ./fargo3d -p -S 20 -o "$opts" "setups/$setup/$setup.par" \
    2>&1 | tee -a "$out/dust_run.log"
else
  echo "=== $CASE St=$ST: restart from output $restart, job $SLURM_JOB_ID, $(date --iso-8601=seconds)" | tee -a "$out/launch.log"
  mpirun -np "$SLURM_NTASKS" ./fargo3d -S "$restart" -o "$opts" "setups/$setup/$setup.par" \
    2>&1 | tee -a "$out/dust_run.log"
fi

if [[ -f "$out/summary$final.dat" ]]; then
  echo "=== $CASE St=$ST complete, $(date --iso-8601=seconds)" | tee -a "$out/launch.log"
else
  echo "=== segment ended at output $(last_complete_output "$out" 20 "$out/gasdens20.dat"), $(date --iso-8601=seconds)" | tee -a "$out/launch.log"
fi
