"""Bed position and vertical domain dimensions for the four simulated slopes.

Adapted from the thesis of F. N. Schuch (04-figuras/Python/Domain-config.py),
self-contained version. Regenerate the figure with:

    uv run --with numpy --with matplotlib python docs/methodology/figures/domain_dimensions.py
"""

import numpy as np
import matplotlib.pyplot as plt

# Number of mesh points supported by Incompact3d (from Domain-config.py)
POINTS = [
    3, 5, 7, 9, 11, 13, 17, 19, 21, 25, 31, 33, 37, 41, 49, 51, 55, 61, 65,
    73, 81, 91, 97, 101, 109, 121, 129, 145, 151, 161, 163, 181, 193, 201,
    217, 241, 251, 257, 271, 289, 301, 321, 325, 361, 385, 401, 433, 451,
    481, 487, 501, 513, 541, 577, 601, 641, 649, 721, 751, 769, 801,
]

XLX = 337.5  # streamwise domain length L_1a + L_1 + L_1b
X0 = 25.0  # inflow sponge zone length L_1a
NX = 1081  # streamwise mesh points n_1
DY = 0.140625  # vertical mesh spacing Delta x_2
T = 297.5  # position where the bed becomes flat, upstream of the outflow
SLOPES = [1.25, 2.5, 5.0, 10.0]  # bed slope S [%]


def domain(s, xlx=XLX, x0=X0, nx=NX, dy=DY, t=T):
    """Bed position and vertical dimensions for a given bed slope s [%]."""
    x = np.linspace(-x0, xlx - x0, num=nx, endpoint=True)
    # minimum domain height: inflow depth + drop over the ramp + bottom layer
    yly_min = 1.0 + s * t / 100.0 + 1.0
    ny = next(n for n in POINTS if (n - 1) * dy >= yly_min)
    yly = dy * float(ny - 1)
    layer = yly - (yly_min - 1.0)  # bottom layer thickness L_2b
    ramp = np.clip(1.0 + s * x / 100.0, 1.0, yly - layer)
    return x, ramp, yly, layer, ny


def main():
    fig, ax = plt.subplots(figsize=(6.0, 3.5), layout="constrained")

    yly_ticks, colors = [], []
    for i, s in enumerate(SLOPES):
        x, ramp, yly, _, _ = domain(s)
        ax.plot(x, ramp, label=rf"$S={s}\%$")
        yly_ticks.append(yly)
        colors.append(f"C{i}")

    ax.set_xlabel(r"$x_1$")
    ax.set_ylabel(r"$x_2$")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(35.0, 0.0)
    ax.axvline(T, linestyle="--", color="silver")
    ax.legend(title="Bed position")

    ax2 = ax.twinx()
    ax2.set_ylim(35.0, 0.0)
    ax2.set_yticks(yly_ticks)
    ax2.set_ylabel(r"$L_2 + L_{2b}$")
    ax2.grid(False)
    for ticklabel, tickcolor in zip(ax2.get_yticklabels(), colors):
        ticklabel.set_color(tickcolor)

    fig.savefig("domain-dimensions.svg")

    print("| $S$ [%] | $L_2$ | $L_{2b}$ | $n_2$ |")
    print("| ------- | ----- | -------- | ----- |")
    for s in SLOPES:
        _, _, yly, layer, ny = domain(s)
        print(f"| {s} | {yly - layer:.6g} | {layer:.6g} | {ny} |")


if __name__ == "__main__":
    main()
