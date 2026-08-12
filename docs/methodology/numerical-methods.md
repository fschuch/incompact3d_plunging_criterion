# Numerical Methodology

The numerical simulations were carried out by
[Incompact3d](https://github.com/xcompact3d/Incompact3d), an open
source tool based on a Boussinesq system for incompressible fluids,
designed for supercomputers. A forked version of the original code was
developed for this work {cite:p}`schuch2020thesis` in order to
simulate the plunging flow, and it is included in this repository (see
[Running the Solver](../running/solver.md)).

The code is based on high-order finite-difference schemes and a
predictor-corrector method for the Navier-Stokes equations, which
leads to the resolution of a Poisson equation for the pressure. With a
spectral approach using three-dimensional (3D) Fast Fourier Transforms
(FFTs) and the concept of modified wavenumbers
{cite:p}`lele1992compact` for this equation, the divergence-free
condition is ensured up to machine accuracy. In order to prevent
spurious pressure oscillations, the pressure mesh is staggered from
the other variables by half a mesh. More information regarding the
numerical aspects of Incompact3d is available in
{cite:t}`laizet2009high`.

The computational domain is divided into several pencils, each one
solved by a different MPI process, resulting in a highly scalable 2D
domain decomposition. A complete description concerning the parallel
strategy can be found in {cite:t}`laizet2011incompact3d`.

The governing equations {eq}`eq-continuity`-{eq}`eq-scalar-transport`
are solved under an implicit Large-Eddy Simulation (iLES) approach, in
which only the largest and energy-containing flow structures are
resolved. The small scales, on the other hand, are not resolved;
instead, they are modelled via an artificial dissipation that takes
place when computing the viscous term
{cite:p}`sagaut2006large,grinstein2007implicit`.
{cite:t}`lamballais2011straightforward` and
{cite:t}`dairay2017numerical` presented a new set of coefficients for
the sixth-order compact finite-difference schemes, designed to be
over-dissipative at the highest wavenumbers and control the aliasing
errors via the viscous term. Notice that the solution only differs
from the classic high-order finite-difference scheme at the highest
wavenumbers, where even the latter becomes inaccurate.

## Simulated cases

Twelve simulations were conducted, aiming to validate the predicted
plunge depth:

| Case    | Bed slope [%] | Settling velocity [-] | Initial densimetric Froude number [-] |
| ------- | ------------- | --------------------- | ------------------------------------- |
| 1.25-0  | 1.25          | 0.0                   | 2.19                                  |
| 1.25-15 | 1.25          | 0.0015                | 1.96                                  |
| 1.25-30 | 1.25          | 0.0030                | 1.75                                  |
| 2.5-0   | 2.50          | 0.0                   | 4.66                                  |
| 2.5-15  | 2.50          | 0.0015                | 4.16                                  |
| 2.5-30  | 2.50          | 0.0030                | 3.72                                  |
| 5.0-0   | 5.00          | 0.0                   | 11.15                                 |
| 5.0-15  | 5.00          | 0.0015                | 9.97                                  |
| 5.0-30  | 5.00          | 0.0030                | 8.90                                  |
| 10.0-0  | 10.00         | 0.0                   | 28.80                                 |
| 10.0-15 | 10.00         | 0.0015                | 25.74                                 |
| 10.0-30 | 10.00         | 0.0030                | 23.00                                 |

The output of all cases is available on Zenodo, see the
[dataset description](../dataset/description.md).
