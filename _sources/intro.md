# Plunging condition for particle-laden flows over sloping bottoms

This website is the companion material for the article
**"Plunging condition for particle-laden flows over sloping bottoms:
three-dimensional turbulence-resolving simulations"**
(Computers & Geosciences, 2021) {cite:p}`schuch2021plunging`, by:

> **Felipe Nornberg Schuch**
> *School of Technology, Pontifical Catholic University of Rio Grande do Sul, Porto Alegre, Brazil.*
>
> **Eckart Meiburg**
> *Department of Mechanical Engineering, University of California Santa Barbara, Santa Barbara, USA.*
>
> **Jorge Hugo Silvestrini**
> *School of Technology, Pontifical Catholic University of Rio Grande do Sul, Porto Alegre, Brazil.*

**Abstract:** Hyperpycnal flows are observed when the density of a fluid entering into a quiescent basin is greater than that of the ambient fluid. This difference can be due to temperature, salinity, turbidity, concentration, or a combination of them. Over a sloping bottom, the inflowing momentum decreases progressively until a critical point is reached where the inflow plunges under the ambient fluid and flows along the bed as an underflow density current. In the present work, a new equation is proposed in order to predict the critical depth for plunging, i.e., the plunging criterion. It differs from previous studies since it includes the role of the settling velocity and the bed slope. The high spatiotemporal resolution from twelve original numerical simulations allows us to validate the initial hypotheses established, in addition to numerical and experimental data available in the literature, and good agreement is found between them. A negative value for the mixing coefficient was observed for the first time for the hyperpycnal flow in a tilted channel. This indicates that if the settling velocity of the suspended material is high enough, the submerged flow may lose fluid to the environment (detrainment), instead of incorporating it. The proposed plunging criterion may assist in the design of future experimental or numerical works.

**Keywords:** Plunging flow, Plunging criterion, Turbidity current, Large-Eddy Simulation.

## What you will find here

::::{grid} 1 1 2 2
:class-container: text-center
:gutter: 3

:::{grid-item-card}
:link: methodology/governing-equations
:link-type: doc
:class-header: bg-light

Methodology 🌊
^^^

The mathematical description of the problem: governing equations,
the computational domain and its boundary conditions, and the
numerical methods employed by the solver.

:::

:::{grid-item-card}
:link: dataset/description
:link-type: doc
:class-header: bg-light

Dataset 📊
^^^

The output from the twelve simulated cases is available on Zenodo in
NetCDF format, together with example notebooks showing how to load,
process and plot it.

:::

:::{grid-item-card}
:link: running/solver
:link-type: doc
:class-header: bg-light

Navier-Stokes Solver 💻
^^^

The simulations were carried out by a customized version of the
open source solver Incompact3d, included in this repository.

:::

:::{grid-item-card}
:link: running/containers
:link-type: doc
:class-header: bg-light

Containers 📦
^^^

Dev container and Apptainer recipes to build and run the Fortran
solver without installing its toolchain on the host system.

:::

::::

```{note}
The results and their discussion are presented in the article itself
and are not reproduced here for copyright reasons. See the
[published version](https://doi.org/10.1016/j.cageo.2021.104880).
The methodology is described in further detail in the PhD thesis of
{cite:t}`schuch2020thesis`.
```

## How to cite this work

```bibtex
@article{schuch2021plunging,
  title = {Plunging criterion for particle-laden flows over sloping bottoms:
           Three-dimensional turbulence-resolving simulations},
  author = {Schuch, Felipe N. and Meiburg, Eckart and Silvestrini, Jorge H.},
  journal = {Computers \& Geosciences},
  volume = {156},
  pages = {104880},
  year = {2021},
  publisher = {Elsevier},
  doi = {10.1016/j.cageo.2021.104880}
}
```
