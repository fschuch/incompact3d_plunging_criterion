# Navier-Stokes Solver

The numerical simulations were carried out by
[Incompact3d](https://github.com/xcompact3d/Incompact3d), an open
source tool based on a Boussinesq system for incompressible fluids,
designed for supercomputers. A forked version of the original code was
developed for this work in order to simulate the plunging flow, and it
is included in the
[`incompact3d`](https://github.com/fschuch/incompact3d_plunging_criterion/tree/main/incompact3d)
folder of this repository. The numerical schemes it employs are
described in the [numerical methodology](../methodology/numerical-methods.md).

## Building natively

The solver requires a Fortran compiler (e.g. `gfortran`), an MPI
implementation (e.g. OpenMPI) and `make`. From the repository root:

```bash
make solver     # same as: make -C incompact3d
```

See `incompact3d/Makefile` for the available build options. If you
prefer not to install the toolchain on your host system, use one of
the provided [containers](containers.md).

## Running a case

The simulation parameters are defined in a `.prm` file (see
`incompact3d/cases/`). To run a case with MPI:

```bash
mpirun -np 4 incompact3d -i BC-Plumes.prm
```

A quick 5-iteration verification case on 2 MPI ranks is available via:

```bash
make smoke-test
```

## License

Incompact3d is distributed under the
[GNU General Public License v3.0](https://github.com/fschuch/incompact3d_plunging_criterion/blob/main/incompact3d/LICENSE).
