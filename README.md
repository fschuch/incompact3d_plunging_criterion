# Plunging condition for particle-laden flows over sloping bottoms: three-dimensional turbulence-resolving simulations

This repository contains the source code employed for our work **"Plunging condition for particle-laden flows over sloping bottoms: three-dimensional turbulence-resolving simulations"** (Computers & Geosciences, 2021), by:

> **Felipe Nornberg Schuch**  
> *School of Technology, Pontifical Catholic University of Rio Grande do Sul, Porto Alegre, Brazil.*  
> felipe.schuch@edu.pucrs.br
>  
> **Eckart Meiburg**  
> *Department of Mechanical Engineering, University of California Santa Barbara, Santa Barbara, USA*  
> meiburg@engineering.ucsb.edu 
>
> **Jorge Hugo Silvestrini**  
> *School of Technology, Pontifical Catholic University of Rio Grande do Sul, Porto Alegre, Brazil.*  
> jorgehs@pucrs.br
>
> **Abstract:** Hyperpycnal flows are observed when the density of a fluid entering into a quiescent basin is greater than that of the ambient fluid. This difference can be due to temperature, salinity, turbidity, concentration, or a combination of them. Over a sloping bottom, the inflowing momentum decreases progressively until a critical point is reached where the inflow plunges under the ambient fluid and flows along the bed as an underflow density current. In the present work, a new equation is proposed in order to predict the critical depth for plunging, i.e., the plunging criterion. It differs from previous studies since it includes the role of the settling velocity and the bed slope. The high spatiotemporal resolution from twelve original numerical simulations allows us to validate the initial hypotheses established, in addition to numerical and experimental data available in the literature, and good agreement is found between them. A negative value for the mixing coefficient was observed for the first time for the hyperpycnal flow in a tilted channel. This indicates that if the settling velocity of the suspended material is high enough, the submerged flow may lose fluid to the environment (detrainment), instead of incorporating it. The proposed plunging criterion may assist in the design of future experimental or numerical works.
>
> **Keywords:** Plunging flow, Plunging criterion, Turbidity current, Large-Eddy Simulation.

## Navier-Stokes Solver

The numerical simulations were carried out by [incompact3d](./incompact3d), an open source tool based on Boussinesq system for incompressible fluids, designed for supercomputers. A forked version from the original code was developed for this work in order to simulate the plunging flow.

## Simulations

Twelve simulations were conducted, aiming to validate the predicted plunge depth:

| Cases | Bed slope [%] | Settling velocity [-] | Initial densimetric Frounde number [-] |
|-------|---------------|-----------------------|----------------------------------------|
| 1.25-0 |   1.25 | 0.0    |    2.19 |
| 1.25-15 |  1.25 | 0.0015 |    1.96 |
| 1.25-30 |  1.25 | 0.0030 |    1.75 |
| 2.5-0 |    2.50 | 0.0    |    4.66 |
| 2.5-15 |   2.50 | 0.0015 |    4.16 |
| 2.5-30 |   2.50 | 0.0030 |    3.72 |
| 5.0-0 |    5.00 | 0.0    |   11.15 |
| 5.0-15 |   5.00 | 0.0015 |    9.97 |
| 5.0-30 |   5.00 | 0.0030 |    8.90 |
| 10.0-0 |  10.00 | 0.0    |   28.80 |
| 10.0-15 | 10.00 | 0.0015 |   25.74 |
| 10.0-30 | 10.00 | 0.0030 |   23.00 |

## Dataset [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4044388.svg)](https://doi.org/10.5281/zenodo.4044388)

Data from the Twelve simulations are included. The output files from incompact3d were converted to [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) aiming to be more friendly than raw binaries, they include data arrays together with metadata and the coordinates:

* `t`: Time;
* `x`: Streamwise coordinate;
* `y`: Vertical coordinate;
* `z`: Spanwise coordinate;

Due to limitations at the storage size, the time sampling is smaller for larger arrays. There are four files for each simulated case:

* `3d-case-<num>.nc`: Complete 3D snapshot at 12 different dimensionless times (250, 500, 750, 1,000, 1,250, 1,500, 1,750, 2,000, 3,000, 4,000, 5,000, 6,000), the arrays are:
  - `ux (x,y,z,t)`: Streamwise velocity;
  - `uy (x,y,z,t)`: Vertical velocity;
  - `uz (x,y,z,t)`: Spanwise velocity;
  - `phi (x,y,z,t)`: Concentration.
* `xy-planes-case-<num>.nc`: Spanwise-averaged quantities with a timestep of 5 dimensionless units:
  - `ux (x,y,t)`: Streamwise velocity;
  - `uy (x,y,t)`: Vertical velocity;
  - `uz (x,y,t)`: Spanwise velocity;
  - `phi (x,y,t)`: Concentration.
* `LA-case-<num>.nc`: The complete spatio-temporal analysis of the relevant quantities is possible in a layer-averaged context per width unit, it is based on three variables:
  - `Layer-averaged Uh (x,t)`;
  - `Layer-averaged U2h (x,t)`;
  - `Layer-averaged UCh (x,t)`;

  that can be used to compute layer-averaged velocity `U = U2h/Uh`, flow depth `H = (Uh)²/U2h`, flow discharge `Q = Uh`, concentration `C = UCh/Uh` and local densimetric Froude number `Fr`. In addition to:
  - `utau (x,t)`: Spanwise-averaged bed shear velocity;
  - `dep (x,t)`: Spanwise-averaged deposition rate.

  All arrays with a timestep of 2.5 dimensionless units.

* `steady-state-case-<num>.nc`: Time (last 2,000 dimensionless units for all cases) and spanwise averaged quantities:
  - `ux (x,y)`: Streamwise velocity;
  - `uy (x,y)`: Vertical velocity;
  - `phi (x,y)`: Concentration.

  Together with layer-averaged quantities (`Uh`, `U2h` and `UCh`), bed shear velocity `utau` and deposition rate `dep`, as describe previously.

Each file can be loaded with the Python package [xarray](http://xarray.pydata.org/en/stable/) (see [Why xarray](http://xarray.pydata.org/en/stable/why-xarray.html)), for instance:

```python
dataset = xr.load_dataset("<filename>")
```

## Examples

* [00-LA-and-Convert-to-NetCDF.ipynb](http://nbviewer.jupyter.org/github/fschuch/incompact3d_plunging_criterion/blob/main/Notebooks/00-LA-and-Convert-to-NetCDF.ipynb) - This Notebook is presented only for reference, it was used to compute the layer-averaged quantities and to convert the raw binary data from Xcompact3d to NetCDF. These files are available at [Zenodo](https://doi.org/10.5281/zenodo.3968993);

* [01-Computing-and-Plotting.ipynb](http://nbviewer.jupyter.org/github/fschuch/incompact3d_plunging_criterion/blob/main/Notebooks/01-Computing-and-Plotting.ipynb) - This Notebook shows how to read, compute and plot the variables presented in our work.

## Setup

The examples above are only rendered online, if you choose to install and run them locally, follow these steps:

1. Clone the repository:
```bash
git clone https://github.com/fschuch/incompact3d_plunging_criterion.git
```

2. Install the environment. This repository includes an [`environment.yaml`](environment.yaml) file containing a list of the requirements to run the examples. To install them with [`Anaconda`](https://www.anaconda.com/) run:
```bash
conda env create -f environment.yml
conda activate plunging-criterion
```

3. Download the data files from [Zenodo](https://doi.org/10.5281/zenodo.4044388) to the root folder of this repository.

4. Start a Jupyter session:
```bash
jupyter lab
```

## Copyright and License

Incompact3d is distributed under [GNU General Public License v3.0](./incompact3d/LICENSE).

Copyright (c) 2021, Felipe N. Schuch. The remaining content in this repository is under Creative Commons Attribution [CC-BY 4.0](https://opensource.org/licenses/BSD-3-Clause).



