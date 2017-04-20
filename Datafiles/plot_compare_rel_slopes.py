#!/usr/bin/env python3
import itertools, sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import utils

with utils.plot(__file__):
    d = utils.load_table("compare_rel_slopes.txt", na_filter=False)

    # we want to use (num_particles, freq) as an axis
    # but matplotlib doesn't let us do that
    # so we have to use this awful hack
    interactions = [
        ("", r"$(\sigma_{\mathrm{A}} = 0, \sigma_{\mathrm{B}} = \infty)$"),
        ("_sigmaB=4.0", r"$(\sigma_{\mathrm{A}} = 0, \sigma_{\mathrm{B}} = 4)$"),
        ("_sigmaA=0.5", r"$(\sigma_{\mathrm{A}} = \frac{1}{2}, \sigma_{\mathrm{B}} = \infty)$"),
        ("_sigmaA=0.5_sigmaB=4.0", r"$(\sigma_{\mathrm{A}} = \frac{1}{2}, \sigma_{\mathrm{B}} = 4)$"),
    ]
    interaction_colors = {
        "": "#a883e4",
        "_sigmaB=4.0": "#951c16",
        "_sigmaA=0.5": "#005e93",
        "_sigmaA=0.5_sigmaB=4.0": "#1e1e1e",
    }
    interaction_markers = {
        "": "o",
        "_sigmaB=4.0": "X",
        "_sigmaA=0.5": "o",
        "_sigmaA=0.5_sigmaB=4.0": "x",
    }
    interaction_markerfacecolors = {
        "": None,
        "_sigmaB=4.0": None,
        "_sigmaA=0.5": "none",
        "_sigmaA=0.5_sigmaB=4.0": "none",
    }
    interaction_label = dict(interactions)
    interaction_pos = dict((x, i) for i, (x, _) in enumerate(interactions))
    systems = [(num_particles, freq) for freq, num_particles in
               itertools.product([1.0, 0.28], [6, 12, 20])]
    system_pos = dict((x, i) for i, x in enumerate(systems))
    d["system_pos"] = d.apply(
        lambda r: system_pos[(r["num_particles"], r["freq"])], axis=1)

    fig, ax = plt.subplots(figsize=(5, 3.5))
    for inter, g in sorted(d.groupby("interaction"),
                           key=lambda x: interaction_pos[x[0]]):
        g.sort_values(["freq", "num_particles"], inplace=True,
                      ascending=[False, True])
        ax.plot(g["system_pos"],
                abs(g["slope"]),
                linestyle="",
                linewidth=2,
                marker=interaction_markers[inter],
                markerfacecolor=interaction_markerfacecolors[inter],
                color=interaction_colors[inter],
                label=interaction_label[inter])

    ax.set_xticks(range(len(systems)))
    ax.set_xticklabels(systems, rotation=90)
    ax.set_xlabel(r"$(N, \omega)$")
    ax.set_ylabel(r"$|\rho_{15}|$")
    ax.set_yscale("log")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    fig.subplots_adjust(right=0.65)
    utils.savefig(fig, fn="../Manuscript/fig-rel-slopes2.pdf")
