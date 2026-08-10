FROM ubuntu:24.04

# Toolchain for the legacy Incompact3d Fortran code (Makefile defaults:
# i3d_CMP=gcc -> mpif90, i3d_FFT=generic). FFTW is included so the
# fftw3/fftw3_f03 FFT backends can also be built if needed.
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        make \
        git \
        openmpi-bin \
        libopenmpi-dev \
        libfftw3-dev \
        gdb \
        less \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Allow mpirun as root inside the container (handy when the container user is root)
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
