# FARGO3D for Yu, Li & Lee (2026)

The FARGO3D code of Yu, Li & Lee (2026). It is public FARGO3D (`f359328`,
December 2024) with the additions listed below and the setups of the paper.

FARGO3D documentation: <https://fargo3d.github.io/documentation>.

## Additions to FARGO3D

- **Dust diffusivity.** `D = alpha c_s^2 / Omega_K / (1 + St^2)`, including
  the colatitude faces of an adiabatic disk, with an optional reduction by
  `1 + rho_d/rho_g` (`DUSTDIFFD2G_REDUCTION`).
- **Beta cooling** (`BETACOOLING`). The gas relaxes toward its initial
  temperature on `t_cool = Beta / Omega_K`.
- **Damping zones.** At the inner and outer radial edges, the gas is damped
  toward its initial state and the dust toward its steady drift solution.
- **Setups of the paper.** Initial state, analytic boundaries (in
  `std/boundaries.txt`), gas and dust floors, the dust velocity caps of the
  Athena++ problem generator (`setups/modelB_*/floor.cu`, GPU only), and the
  dust release.
- **Build.** Python 3.13 and CUDA 12 support, `FARGO_CUDA_ARCH` passed to
  nvcc, one GPU per local MPI rank, and a setup's handwritten `floor.cu`
  always compiled as written.

## Setups

The disk follows Huang et al. (2025): `h = 0.05`, `alpha = 1e-2` and
`Beta = 1`, with a planet of 2.4 thermal masses (mass ratio `3e-4`) on a fixed
circular orbit at `R0 = 1`, on a 768 x 256 x 128 spherical mesh. A run has two
stages:

| setup | stage |
|---|---|
| `modelB_gas` | gas only, 0 to 100 orbits (outputs 0 to 20) |
| `modelB_contdiff` | dust released at output 20, with dust diffusion (ContDiff), to 1500 orbits |
| `modelB_nodiff` | the same without dust diffusion (NoDiff) |

The gas stage is the same for every run. The dust stages take the Stokes
number from `Invstokes1`: 10, 31.62 and 100 for St = 0.1, 0.0316 and 0.01
(the parameter files ship with 10).

## Running

Build each setup, then copy the executable next to its parameter file:

```bash
make SETUP=modelB_gas UNITS=0 RESCALE=0 GPU=1 PARALLEL=1 MPICUDA=1
cp fargo3d setups/modelB_gas/
make mrproper
make SETUP=modelB_contdiff UNITS=0 RESCALE=0 GPU=1 PARALLEL=1 MPICUDA=1
cp fargo3d setups/modelB_contdiff/
```

Run the gas stage once:

```bash
mpirun -np 4 setups/modelB_gas/fargo3d setups/modelB_gas/modelB_gas.par
```

Each dust run then starts from a copy of its output 20. With `-p` the dust is
released, e.g. for St = 0.0316 with diffusion:

```bash
mkdir -p outputs/contdiff_St0.0316
cp outputs/modelB_gas/*20.dat outputs/modelB_gas/planet0.dat outputs/modelB_gas/domain_*.dat outputs/contdiff_St0.0316/
mpirun -np 4 setups/modelB_contdiff/fargo3d -p -S 20 \
  -o "Invstokes1,31.62,OutputDir,@outputs/contdiff_St0.0316" \
  setups/modelB_contdiff/modelB_contdiff.par
```

A dust stage takes about 60 hours on four H100 GPUs. Leave `FARGO_CUDA_ARCH`
unset to build for nvcc's default target, as was done for the paper.

------------------------

Report bugs in FARGO3D itself to the [FARGO3D issues
page](https://github.com/FARGO3D/fargo3d/issues).
