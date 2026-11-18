#!/bin/bash
# Shared functions of the Model B launchers. Sourced by gas_body.sh and dust_body.sh.

set -euo pipefail

FIELDS=(gasdens gasenergy gasvx gasvy gasvz dust1dens dust1vx dust1vy dust1vz)

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
output_root=${FARGO_OUTPUT_ROOT:-"$SCRATCH/fargo3d_rerun"}

: "${SLURM_JOB_ID:?Run this inside a Slurm job.}"
: "${SLURM_TMPDIR:?SLURM_TMPDIR is not set.}"
if [[ ${SLURM_NTASKS:-0} -ne 4 ]]; then
  echo "ERROR: expected four MPI ranks (one per H100), got SLURM_NTASKS=${SLURM_NTASKS:-unset}." >&2
  exit 2
fi

load_modules() {
  module --force purge
  # shellcheck disable=SC2086
  module load $FARGO_MODULES
  command -v nvcc mpirun python3
  export OMP_NUM_THREADS=1
  export OMPI_MCA_opal_cuda_support=true
  export FARGO_CUDA_ARCH
}

# Build a setup in a private copy of the repository, so that concurrent jobs
# never share bin/ or the executable. Prints the build directory.
build_setup() {
  local setup=$1 logdir=$2
  local bdir="$SLURM_TMPDIR/fargo_$setup"
  rsync -a --exclude .git --exclude outputs "$repo_root/" "$bdir/"
  # FARGO3D's make must run serially: its code generators share bin/.
  if ! (cd "$bdir" && make SETUP="$setup" UNITS=0 RESCALE=0 GPU=1 PARALLEL=1 MPICUDA=1) \
       > "$logdir/build_$setup.log" 2>&1; then
    echo "ERROR: build of $setup failed; see $logdir/build_$setup.log" >&2
    exit 3
  fi
  test -x "$bdir/fargo3d"
  echo "$bdir"
}

# Record where and with what a job segment ran (appended per segment).
write_provenance() {
  local dir=$1 bdir=$2
  {
    echo "=== segment: job $SLURM_JOB_ID on $FARGO_CLUSTER_NAME ($(hostname)), $(date --iso-8601=seconds)"
    echo "git_commit=$(git -C "$repo_root" rev-parse HEAD)"
    echo "git_status:"; git -C "$repo_root" status --short
    echo "fargo_cuda_arch=$FARGO_CUDA_ARCH"
    echo "modules=$FARGO_MODULES"
    echo "binary_sha256=$(sha256sum "$bdir/fargo3d" | cut -d' ' -f1)"
    srun --ntasks=1 nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader
  } >> "$dir/provenance.txt" 2>&1
}

# Highest output number n >= $2 whose nine fields all exist with the size of
# the reference field file $3 and whose summary file exists; empty if none.
last_complete_output() {
  local dir=$1 nmin=$2 ref=$3 size n f ok
  size=$(stat -c %s "$ref")
  for n in $(ls "$dir" 2>/dev/null | sed -n 's/^summary\([0-9][0-9]*\)\.dat$/\1/p' | sort -rn); do
    (( n >= nmin )) || continue
    ok=1
    for f in "${FIELDS[@]}"; do
      [[ -f "$dir/$f$n.dat" && $(stat -c %s "$dir/$f$n.dat") -eq $size ]] || { ok=0; break; }
    done
    if (( ok )); then echo "$n"; return 0; fi
  done
  return 0
}
