# Dataset

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4044388.svg)](https://doi.org/10.5281/zenodo.4044388)

Data from the twelve simulations are available on
[Zenodo](https://doi.org/10.5281/zenodo.4044388). The output files
from Incompact3d were converted to
[NetCDF](https://www.unidata.ucar.edu/software/netcdf/), aiming to be
more friendly than raw binaries; they include data arrays together
with metadata and the coordinates:

* `t`: Time;
* `x`: Streamwise coordinate;
* `y`: Vertical coordinate;
* `z`: Spanwise coordinate.

Due to limitations on the storage size, the time sampling is smaller
for larger arrays. There are four files for each simulated case:

* `3d-case-<num>.nc`: Complete 3D snapshot at 12 different
  dimensionless times (250, 500, 750, 1,000, 1,250, 1,500, 1,750,
  2,000, 3,000, 4,000, 5,000, 6,000), the arrays are:
  - `ux (x,y,z,t)`: Streamwise velocity;
  - `uy (x,y,z,t)`: Vertical velocity;
  - `uz (x,y,z,t)`: Spanwise velocity;
  - `phi (x,y,z,t)`: Concentration.
* `xy-planes-case-<num>.nc`: Spanwise-averaged quantities with a
  timestep of 5 dimensionless units:
  - `ux (x,y,t)`: Streamwise velocity;
  - `uy (x,y,t)`: Vertical velocity;
  - `uz (x,y,t)`: Spanwise velocity;
  - `phi (x,y,t)`: Concentration.
* `LA-case-<num>.nc`: The complete spatio-temporal analysis of the
  relevant quantities is possible in a layer-averaged context per
  width unit, it is based on three variables:
  - `Layer-averaged Uh (x,t)`;
  - `Layer-averaged U2h (x,t)`;
  - `Layer-averaged UCh (x,t)`;

  that can be used to compute layer-averaged velocity `U = U2h/Uh`,
  flow depth `H = (Uh)²/U2h`, flow discharge `Q = Uh`, concentration
  `C = UCh/Uh` and local densimetric Froude number `Fr`. In addition
  to:
  - `utau (x,t)`: Spanwise-averaged bed shear velocity;
  - `dep (x,t)`: Spanwise-averaged deposition rate.

  All arrays with a timestep of 2.5 dimensionless units.

* `steady-state-case-<num>.nc`: Time (last 2,000 dimensionless units
  for all cases) and spanwise averaged quantities:
  - `ux (x,y)`: Streamwise velocity;
  - `uy (x,y)`: Vertical velocity;
  - `phi (x,y)`: Concentration.

  Together with layer-averaged quantities (`Uh`, `U2h` and `UCh`),
  bed shear velocity `utau` and deposition rate `dep`, as described
  previously.

Each file can be loaded with the Python package
[xarray](https://xarray.pydata.org/en/stable/) (see
[Why xarray](https://xarray.pydata.org/en/stable/why-xarray.html)),
for instance:

```python
import xarray as xr

dataset = xr.load_dataset("<filename>")
```

## Example notebooks

Two notebooks exemplify how the dataset was produced and how to work
with it:

* [Layer-Averaged quantities and conversion to NetCDF](../examples/00-LA-and-Convert-to-NetCDF.ipynb) -
  presented only for reference, it was used to compute the
  layer-averaged quantities and to convert the raw binary data from
  Incompact3d to NetCDF;
* [Computing and Plotting](../examples/01-Computing-and-Plotting.ipynb) -
  shows how to read, compute and plot the variables presented in our
  work.

```{important}
The notebooks are rendered here with their stored outputs, since the
dataset is not included in this repository. To run them locally,
follow the setup instructions below.
```

## Running the notebooks locally

1. Clone the repository:

   ```bash
   git clone https://github.com/fschuch/incompact3d_plunging_criterion.git
   ```

2. Install the environment. This repository includes an
   [`environment.yml`](https://github.com/fschuch/incompact3d_plunging_criterion/blob/main/environment.yml)
   file containing a list of the requirements to run the examples. To
   install them with [Anaconda](https://www.anaconda.com/), run:

   ```bash
   conda env create -f environment.yml
   conda activate plunging-criterion
   ```

3. Download the data files from
   [Zenodo](https://doi.org/10.5281/zenodo.4044388) to the root folder
   of this repository.

4. Start a Jupyter session:

   ```bash
   jupyter lab
   ```
