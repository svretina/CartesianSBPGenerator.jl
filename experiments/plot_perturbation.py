"""Plot matched-grid dimensionless matrices from compare_perturbation.jl."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("data", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)

names = ("direct", "analytical", "conventional")
titles = ("Original direct", "Analytical perturbation", "Conventional baseline")
colors = ("#af3937", "#2a687b", "#68614e")
n = 81
for p in (2, 4, 6):
    prefix = [args.data / f"p{p}_n{n}_{name}" for name in names]
    matrices = [np.loadtxt(f"{path}_D.tsv") for path in prefix]
    rows = [np.loadtxt(f"{path}_rows.tsv") for path in prefix]
    vmax = max(np.abs(D).max() for D in matrices)
    norm = SymLogNorm(linthresh=0.025, vmin=-vmax, vmax=vmax)
    exponents = np.arange(-2, np.floor(np.log10(vmax)) + 1, 2)
    decades = 10.0 ** exponents
    ticks = np.r_[-decades[::-1], 0, decades]
    fig = plt.figure(figsize=(13.2, 12), layout="constrained")
    gs = fig.add_gridspec(4, 6, height_ratios=(1, 1, 0.9, 0.9))
    node = np.arange(n)
    heat_axes = []
    for j, (title, D) in enumerate(zip(titles, matrices)):
        nz = fig.add_subplot(gs[0, 2*j:2*j+2])
        nz.spy(D, markersize=1.6, color=colors[j])
        nz.set(title=title, ylabel="Row" if j == 0 else None)
        nz.set_xlabel("Column", labelpad=-2)
        heat = fig.add_subplot(gs[1, 2*j:2*j+2])
        image = heat.imshow(D, norm=norm, cmap="RdBu_r", interpolation="none")
        heat.set(ylabel="Row" if j == 0 else None, xlabel="Column")
        heat_axes.append(heat)
    fig.colorbar(image, ax=heat_axes, shrink=0.9, ticks=ticks,
                 label="Dimensionless coefficient hD")

    for panel, (i, desc) in enumerate(((0, "Endpoint"), (4, "Fifth row"),
                                      (n//2, "Interior phase A"))):
        ax = fig.add_subplot(gs[2, 2*panel:2*panel+2])
        for label, color, D in zip(titles, colors, matrices):
            ids = np.flatnonzero(D[i])
            ax.plot(ids-i, D[i, ids], "o-", color=color, lw=1, ms=3, label=label)
        ax.axhline(0, color="#aaa", lw=0.6)
        ax.axvline(0, color="#aaa", lw=0.6)
        ax.set(title=desc, xlabel="Offset from evaluation node",
               ylabel="hD" if panel == 0 else None, yscale="symlog")
    mass = fig.add_subplot(gs[3, :3])
    leading = fig.add_subplot(gs[3, 3:])
    for label, color, values in zip(titles, colors, rows):
        mass.plot(node, values[:, 0], color=color, lw=1.2, label=label)
        leading.plot(node, values[:, 1], color=color, lw=1.2)
    mass.set(xlabel="Node index", ylabel="Diagonal norm weight Hii/h")
    leading.set(xlabel="Node index", ylabel="Signed leading Taylor coefficient")
    leading.set_yscale("symlog", linthresh=0.025)
    mass.legend(loc="best", fontsize=8)
    fig.suptitle(f"Order {p} alternating SBP · {n} nodes on [0, 1]", fontsize=15)
    fig.savefig(args.output / f"alternating_locality_p{p}.png", dpi=145)
    plt.close(fig)
