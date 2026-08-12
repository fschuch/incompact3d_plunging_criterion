# Containers

Two container options are provided to build and run the Fortran
solver without installing its toolchain (gfortran, OpenMPI, make) on
the host system.

## Dev container

Open the repository in VS Code (or any Dev Containers-compatible tool)
and choose *Reopen in Container*. The
[.devcontainer](https://github.com/fschuch/incompact3d_plunging_criterion/tree/main/.devcontainer)
configuration ships the full toolchain (gfortran, OpenMPI, FFTW) plus
Apptainer itself. Then, from the repository root:

```bash
make solver     # build the Fortran solver (same as: make -C incompact3d)
```

```bash
make sif        # build the Apptainer image (incompact3d.sif, uses sudo)
```

```bash
make smoke-test # run a quick 5-iteration verification case on 2 MPI ranks
```

## Apptainer

The definition file
[apptainer/incompact3d.def](https://github.com/fschuch/incompact3d_plunging_criterion/blob/main/apptainer/incompact3d.def)
produces a self-contained image with the compiled solver and its own
OpenMPI. Build it with `make sif` (works inside the dev container or
on any Linux machine with Apptainer), or pull the image published by
CI (on tagged releases) with:

```bash
apptainer pull oras://ghcr.io/fschuch/incompact3d_plunging_criterion:latest
```

Run on a single node using the MPI inside the container:

```bash
apptainer exec incompact3d.sif mpirun -np 4 incompact3d -i BC-Plumes.prm
```

or in hybrid mode with the host MPI for multi-node jobs (requires a
host OpenMPI ABI-compatible with the one in the image, i.e. same major
version):

```bash
mpirun -np 128 apptainer exec incompact3d.sif incompact3d -i BC-Plumes.prm
```

```{note}
When Apptainer itself runs inside a Docker container (as in the dev
container), add `--unsquash` to `apptainer exec` — a FUSE limitation
prevents executing binaries from the mounted image otherwise.
`make smoke-test` already does this; it is harmless on regular hosts.
```

## Continuous integration

The GitHub Actions workflow
[apptainer.yml](https://github.com/fschuch/incompact3d_plunging_criterion/blob/main/.github/workflows/apptainer.yml)
builds the image, runs a short smoke-test simulation
([cases/smoke-test](https://github.com/fschuch/incompact3d_plunging_criterion/tree/main/incompact3d/cases/smoke-test))
on 2 MPI ranks, and publishes the image to the GitHub Container
Registry on version tags (`v*`). In compliance with the GPL-3.0
license, the image also carries the full solver source code at
`/opt/incompact3d-src`.
