# Convenience targets for the solver and its Apptainer image.
# Apptainer targets work inside the dev container (which ships apptainer)
# or on any Linux machine with Apptainer installed.

MPI_RANKS ?= 2

.PHONY: all solver sif smoke-test docs docs-serve clean

all: solver

# Build the Fortran solver natively (see incompact3d/Makefile for options)
solver:
	$(MAKE) -C incompact3d

# Build the Apptainer image (building from a definition file needs root)
incompact3d.sif: apptainer/incompact3d.def $(wildcard incompact3d/*.f90) incompact3d/Makefile
	sudo apptainer build $@ $<

sif: incompact3d.sif

# Short verification run: 5 time steps on $(MPI_RANKS) MPI ranks, using the
# coarse case in incompact3d/cases/smoke-test. --unsquash works around a FUSE
# limitation when apptainer runs inside a Docker container (e.g. the dev
# container); it is harmless elsewhere.
smoke-test: incompact3d.sif
	rm -rf smoke-test-run
	mkdir -p smoke-test-run/data
	cp incompact3d/cases/smoke-test/BC-Plumes.prm smoke-test-run/
	cd smoke-test-run && apptainer exec --unsquash ../incompact3d.sif \
		mpirun -np $(MPI_RANKS) incompact3d -i BC-Plumes.prm

# Build the documentation into docs/_build/html. Dependencies are managed
# by uv (https://docs.astral.sh/uv/), pinned by docs/uv.lock for
# reproducibility; the environment lives in docs/.venv.
docs:
	uv run --project docs --frozen jupyter-book build docs

# Live-reloading preview of the documentation
docs-serve: docs
	uv run --project docs --frozen sphinx-autobuild docs docs/_build/html --open-browser

clean:
	$(MAKE) -C incompact3d clean
	rm -rf incompact3d.sif smoke-test-run docs/_build
