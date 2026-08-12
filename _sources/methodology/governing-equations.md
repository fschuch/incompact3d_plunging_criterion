# Governing Equations

```{note}
The methodology chapters are adapted from the PhD thesis of
{cite:t}`schuch2020thesis`, which describes the numerical framework
employed in the article {cite:p}`schuch2021plunging`.
```

The present framework is based on the resolution of the incompressible
Navier-Stokes equations. $N$ scalar transport equations under the
Boussinesq approximation are used for sediments with different grain
sizes. The incompressible Navier-Stokes equations are written in their
dimensionless form as

```{math}
:label: eq-continuity
\frac{\partial u_j}{\partial x_j} = 0,
```

```{math}
:label: eq-ns
\frac{\partial u_i}{\partial t}
+ u_j \frac{\partial u_i}{\partial x_j}
= - \frac{\partial p}{\partial x_i}
+ \frac{1}{Re} \frac{\partial^2 u_i}{\partial x_j \partial x_j}
+ e^g_i \frac{1}{{Fr_0}^2} \sum_{\ell=1}^N c_\ell,
```

```{math}
:label: eq-scalar-transport
\frac{\partial c_\ell}{\partial t}
+ \left( u_j + u_{s,\ell}\, e^g_j \right)
\frac{\partial c_\ell}{\partial x_j}
= \frac{1}{Re\, Sc} \frac{\partial^2 c_\ell}{\partial x_j \partial x_j},
\quad \ell = 1, \dots, N,
```

where $u_i$, $p$ and $c$ correspond to velocity, pressure and portion
of each particle size, respectively, together with the coordinate
system $x_i$ (see the [computational domain](computational-domain.md)),
time $t$, and the unit vector pointing in the gravity direction
$e^g = [0, 1, 0]$. The settling velocity for each particle size is
$u_s$ and may be considered constant. It is related to the particle
diameter by the Stokes settling velocity law {cite:p}`julien95`, which
assumes that the dominant flow force on an individual particle is the
Stokes drag.

The Reynolds, densimetric Froude and Schmidt numbers at the inlet and
the Stokes settling velocity are the four dimensionless parameters,
defined as follows:

```{math}
:label: eq-dimensionless-parameters
\begin{align}
Re &= \dfrac{\tilde{Q}_0}{\tilde{\nu}}, \\
Fr &= \dfrac{\tilde{Q}_0}{\sqrt{R \tilde{C}_0 \tilde{g} \tilde{h}_0^3}}, \\
Sc &= \dfrac{\tilde{\nu}}{\tilde{D}}, \\
u_{s,\ell} &= \dfrac{\tilde{u}_{s,\ell} \tilde{h}_0}{\tilde{Q}_0}
            = \dfrac{\tilde{d}_{s,\ell}^2 R \tilde{g} \tilde{h}_0}
                    {18 \tilde{\nu} \tilde{Q}_0},
\quad \ell = 1, \dots, N,
\end{align}
```

where $\tilde{\cdot}$ corresponds to dimensional quantities.
$\tilde{Q}_0$ is the volumetric discharge per unit width, defined as
the product of the inlet velocity $\tilde{U}_0$ with the initial depth
$\tilde{h}_0$. The kinematic viscosity is $\tilde{\nu}$; the fresh
water and sediment densities are $\tilde{\rho}_w$ and
$\tilde{\rho}_s$, with $R$ defined as
$(\tilde{\rho}_s - \tilde{\rho}_w)/\tilde{\rho}_w$. The gravity
acceleration is $\tilde{g}$, the constant diffusivity of particle
concentration is $\tilde{D}$ and the grain size is $\tilde{d}_\ell$.

All parameters and variables are made dimensionless using the influx
velocity $\tilde{U}_0$, depth $\tilde{h}_0$ and sediment concentration
$\tilde{C}_0$ (the associated dimensionless parameters are $h_0$,
$U_0$ and $C_0$).
