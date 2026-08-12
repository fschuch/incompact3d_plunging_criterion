# Plunging condition for particle-laden flows over sloping bottoms: three-dimensional turbulence-resolving simulations

[![Docs](https://img.shields.io/badge/docs-fschuch.github.io-blue)](https://fschuch.github.io/incompact3d_plunging_criterion)
[![DOI](https://zenodo.org/badge/DOI/10.1016/j.cageo.2021.104880.svg)](https://doi.org/10.1016/j.cageo.2021.104880)
[![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4044388.svg)](https://doi.org/10.5281/zenodo.4044388)

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
> **Abstract:** Hyperpycnal flows are observed when the density of a fluid entering into a quiescent basin is greater than that of the ambient fluid. This difference can be due to temperature, salinity, turbidity, concentration, or a combination of them. Over a sloping bottom, the inflowing momentum decreases progressively until a critical point is reached where the inflow plunges under the ambient fluid and flows along the bed as an underflow density current. In the present work, a new equation is proposed in order to predict the critical depth for plunging, i.e., the plunging criterion. It differs from previous studies since it includes the role of the settling velocity and the bed slope. The high spatiotemporal resolution from twelve original numerical simulations allows us to validate the initial hypotheses established, in addition to numerical and experimental data available in the literature, and good agreement is found between them. A negative value for the mixing coefficient was observed for the first time for the hyperpycnal flow in a tilted channel. This indicates that if the settling velocity of the suspended material is high enough, the submerged flow may lose fluid to the environment (detrainment), instead of incorporating it. The proposed plunging criterion may assist in the design of future experimental or numerical works.
>
> **Keywords:** Plunging flow, Plunging criterion, Turbidity current, Large-Eddy Simulation.

## Documentation

The full companion material is available at
**[fschuch.github.io/incompact3d_plunging_criterion](https://fschuch.github.io/incompact3d_plunging_criterion)**,
including:

* **[Methodology](https://fschuch.github.io/incompact3d_plunging_criterion/methodology/governing-equations.html)** — the governing equations, the mathematical description of the computational domain and its boundary conditions, and the numerical methods;
* **[Dataset](https://fschuch.github.io/incompact3d_plunging_criterion/dataset/description.html)** — the output from the twelve simulated cases, available on [Zenodo](https://doi.org/10.5281/zenodo.4044388) in NetCDF format, together with example notebooks ([docs/examples](./docs/examples)) showing how to load, process and plot it;
* **[Running the solver](https://fschuch.github.io/incompact3d_plunging_criterion/running/solver.html)** — how to build and run the customized [incompact3d](./incompact3d) solver, natively or with the provided [containers](https://fschuch.github.io/incompact3d_plunging_criterion/running/containers.html) (dev container and Apptainer).

## Copyright and License

Incompact3d is distributed under [GNU General Public License v3.0](./incompact3d/LICENSE).

Copyright (c) 2021, Felipe N. Schuch. The remaining content in this repository is under Creative Commons Attribution [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).
