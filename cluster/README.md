# Slurm launchers for the Alliance clusters

These are launchers for Fir, Rorqual, Nibi and Trillium (GPU). Each job:
- uses four H100 GPUs, with one MPI rank per GPU;
- builds its setup in `$SLURM_TMPDIR` with `sm_90`, so concurrent jobs never share a build;
- writes to `$FARGO_OUTPUT_ROOT`, which defaults to `$SCRATCH/fargo3d_rerun`;
- appends the commit, modules and binary hash to `provenance.txt` in the output directory.

Submit from the repository root. On Trillium, submit from the GPU login node `trig-login01`.

| stage | script | output | time |
|---|---|---|---|
| gas, 0-100 orbits (outputs 0-20) | `cluster/<cluster>/gas.sh` | `modelB_gas` | about 4 h |
| dust, 100-1500 orbits (outputs 20-300) | `cluster/<cluster>/dust.sh` with `CASE`, `ST` | `<CASE>_St<ST>` | about 60 h |

`CASE` is `contdiff` or `nodiff`, and `ST` is `0.1`, `0.0316` or `0.01`.

The dust jobs request 24 h. Each job continues from the last complete output, and a job whose run has reached output 300 exits at once. A run therefore takes three 24-hour jobs, chained with `--dependency=afterany`. Shorter jobs usually start much sooner than 72-hour ones (see "Queue estimates" below). On Fir, Rorqual and Nibi, `--time=72:00:00` at submission runs a whole dust stage in one job.

A continuation refuses to start if the repository's commit differs from the one that started the run. Set `ALLOW_COMMIT_CHANGE=1` to override.

## A test first: gas stage and one dust run

```bash
g=$(sbatch --parsable cluster/fir/gas.sh)
d1=$(sbatch --parsable --dependency=afterok:$g --job-name=cd-St0.1 --export=ALL,CASE=contdiff,ST=0.1 cluster/fir/dust.sh)
d2=$(sbatch --parsable --dependency=afterany:$d1 --job-name=cd-St0.1 --export=ALL,CASE=contdiff,ST=0.1 cluster/fir/dust.sh)
sbatch --dependency=afterany:$d2 --job-name=cd-St0.1 --export=ALL,CASE=contdiff,ST=0.1 cluster/fir/dust.sh
```

## All six dust runs on one cluster

Replace `fir` with `rorqual`, `nibi` or `trillium` to run on another cluster.

```bash
g=$(sbatch --parsable cluster/fir/gas.sh)
for run in contdiff:0.1 contdiff:0.0316 contdiff:0.01 nodiff:0.1 nodiff:0.0316 nodiff:0.01; do
  c=${run%%:*}; s=${run##*:}; dep=afterok:$g
  for seg in 1 2 3; do
    j=$(sbatch --parsable --dependency=$dep --job-name=${c:0:2}-St$s --export=ALL,CASE=$c,ST=$s cluster/fir/dust.sh)
    dep=afterany:$j
  done
done
```

To split the runs over clusters, run the gas stage once and copy its output 20 to the other cluster's output root. On that cluster, submit only the dust jobs, without the dependency on `$g`:

```bash
# on the receiving cluster; /scratch/fyyu is $SCRATCH on Fir (check with `echo $SCRATCH` there)
mkdir -p $SCRATCH/fargo3d_rerun/modelB_gas
rsync -av --include='*20.dat' --include='planet0.dat' --include='domain_*.dat' --exclude='*' \
  fyyu@fir.alliancecan.ca:/scratch/fyyu/fargo3d_rerun/modelB_gas/ $SCRATCH/fargo3d_rerun/modelB_gas/
```

## Queue estimates

`sbatch --test-only` returns the scheduler's current estimate of the start time. It does not submit a job.

```bash
# Fir, Rorqual
for t in 08:00:00 24:00:00 72:00:00; do
  sbatch --test-only --account=def-evelee --time=$t --nodes=1 --ntasks-per-node=4 --cpus-per-task=2 \
    --gpus-per-node=h100:4 --mem=160G --wrap=true 2>&1 | tail -1
done
# Nibi
for t in 08:00:00 24:00:00 72:00:00; do
  sbatch --test-only --account=def-evelee --time=$t --nodes=1 --ntasks=4 --cpus-per-task=4 \
    --gpus-per-node=h100:4 --mem=160G --wrap=true 2>&1 | tail -1
done
# Trillium (trig-login01; 24 h at most)
for t in 08:00:00 24:00:00; do
  sbatch --test-only --account=def-evelee --time=$t --nodes=1 --ntasks-per-node=4 \
    --gpus-per-node=4 --wrap=true 2>&1 | tail -1
done
```
