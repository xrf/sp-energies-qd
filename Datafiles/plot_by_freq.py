#!/usr/bin/env python3
# Generates ../FigureFiles/fig-by-freq-*.svg from the data files.
import argparse, itertools, os, sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import utils

def plot(ax, label, interaction, num_shells, num_filled, title=False):
    d = utils.filter_preferred_ml(utils.load_all())
    d = d[(d["method"] != "imsrg[f]+eom[n]") &
          (d["method"] != "magnus_quads+eom") &
          (d["interaction"] == interaction) &
          (d["label"] == label) &
          (d["num_shells"] == num_shells) &
          (d["num_filled"] == num_filled) &
          # filter out higher frequencies because they stretch horizontal axis
          # of plot too much
          (d["freq"] <= 1.0)]
    base_markersize = 5.0
    for method, case in d.groupby("method"):
        case = case.sort_values("freq").drop_duplicates()
        xs = case["freq"]
        ys_ref = xs.map(lambda freq: d[(d["method"] == "hf") &
                                       (d["freq"] == freq)]["energy"])
        ys = case["energy"] / ys_ref
        marker = utils.METHOD_MARKER[method]
        style = {
            "marker": marker,
            "markerfacecolor": "none",
            "markersize": (utils.MARKERSIZE_CORRECTION.get(marker, 1.0) *
                           base_markersize),
            "color": utils.METHOD_COLOR[method],
            "label": utils.METHOD_LABEL[method],
        }
        ax.plot(xs, ys, **style)
        yield method, matplotlib.lines.Line2D([], [], **style)
    if title:
        ml, = d["ml"].unique()
        energy_type = {"ground": "ground state",
                       "add": "addition", "rm": "removal"}[label]
        ax.set_title("{energy_type} energy for {num_particles} particles "
                     r"($M_\ell = {ml}, K = {num_shells}$)"
                     .format(**locals()))
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"${0} / {0}_{{\mathrm{{HF}}}}$"
                  .format(utils.ENERGY_SYMBOL[label]))
    fig.tight_layout()

num_shells = 10
num_filled = 2
num_particles = num_filled * (num_filled + 1)

p = argparse.ArgumentParser()
p.add_argument("-i", "--interactive", action="store_true")
p.add_argument("-t", "--title", action="store_true")
p.add_argument("--interaction", default="normal")
args = p.parse_args()
interaction = args.interaction

plt.rcParams["interactive"] = args.interactive
with utils.plot(__file__):
    width = 6.5
    height = 4.0
    fig = plt.figure(figsize=(width, height))
    gs = matplotlib.gridspec.GridSpec(2, 2)
    entries = {}

    entries.update(plot(fig.add_subplot(gs[0]), "ground", interaction,
                        num_shells, num_filled))
    entries.update(plot(fig.add_subplot(gs[2]), "add", interaction,
                        num_shells, num_filled))
    entries.update(plot(fig.add_subplot(gs[3]), "rm", interaction,
                        num_shells, num_filled))

    ax = fig.add_subplot(gs[1])
    ax.axis("off")
    ax.legend(handles=[v for _, v in sorted(entries.items())], frameon=False,
              loc="center", bbox_to_anchor=(0.5, 0.5))
    fig.tight_layout()
    utils.savefig(fig,
                  fn=("../Manuscript/fig-by-freq-{num_shells}-{num_particles}-"
                      "{interaction}.pdf".format(**locals())))
