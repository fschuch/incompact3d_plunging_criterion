# Dependencies

Incompact3d is OS independent, however, Fortran-90 and MPI compilers  are requeride.
For instance, in an Ubuntu instalation, the following packages should be installed:

    sudo apt install gfortran libopenmpi-dev

# Compilation

There are two ways to set the necessary compilation flags. First, uncoment the following lines at the [Makefile](./incompact3d/Makefile):

```makefile
i3d_FLOW_TYPE = Plumes
i3d_LCL = local#local,lad,sdu
i3d_IVER = 18#16,17
i3d_CMP = gcc#intel,gcc
i3d_FFT = generic#mkl,generic,fftw3_f03,fftw3
```

It can also be done with environment variables, for instance, setting them at the `~/.bashrc` file:
```
export i3d_FLOW_TYPE="Plumes"
export i3d_LCL="local"
export i3d_IVER="18"
export i3d_CMP="gcc"
export i3d_FFT="generic"
```

The last steep is type:

    make

# Running the code

Now you should be able to run incompact like this:

    mpirun -n 4 ./incompact3d
