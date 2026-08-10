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
        sudo \
        software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Apptainer (from the official PPA: has amd64 and arm64 builds), so the
# apptainer/incompact3d.def image can be built and tested inside the dev
# container. Building a .sif needs privileges: the devcontainer runs with
# --privileged (see devcontainer.json) and the ubuntu user gets passwordless
# sudo, mirroring the CI runner (sudo apptainer build ...).
RUN add-apt-repository -y ppa:apptainer/ppa \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends apptainer \
    && rm -rf /var/lib/apt/lists/* \
    && echo "ubuntu ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/ubuntu

# Allow mpirun as root inside the container (handy when the container user is root)
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
