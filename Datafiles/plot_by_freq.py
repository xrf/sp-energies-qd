#!/usr/bin/env python3
# Generates ../FigureFiles/fig-by-freq-*.svg from the data files.
import argparse, itertools, os, sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import utils

p = argparse.ArgumentParser()
p.add_argument("-i", "--interactive", action="store_true")
p.add_argument("-t", "--title", action="store_true")
p.add_argument("-int", "--interaction", default="normal")
p.add_argument("label", metavar="type", help="ground, add, or rm")
kwargs = vars(p.parse_args())
plt.rcParams["interactive"] = kwargs["interactive"]
with utils.plot(__file__):
    d = utils.load_all()
    num_shells = 10
    num_filled = 2
    interaction = kwargs["interaction"]
    label = kwargs["label"]
    ml = utils.label_num_filled_to_ml(label, num_filled)
    d = d[(d["method"] != "imsrg[f]+eom[n]") &
          (d["method"] != "magnus_quads+eom") &
          (d["interaction"] == interaction) &
          (d["label"] == label) &
          (d["num_shells"] == num_shells) &
          (d["num_filled"] == num_filled) &
          # filter out higher frequencies because they stretch the plot too much
          (d["freq"] <= 1.0) &
          (d["ml"] == ml)]
    num_particles = num_filled * (num_filled + 1)
    energy_type = {"ground": "ground state",
                   "add": "addition", "rm": "removal"}[label]
    fig, ax = plt.subplots()
    fig.set_size_inches((4, 3))
    for method, case in d.groupby("method"):
        case = case.sort_values("freq").drop_duplicates()
        xs = case["freq"]
        ys_ref = xs.map(lambda freq: d[(d["method"] == "hf") &
                                       (d["freq"] == freq)]["energy"])
        ys = case["energy"] / ys_ref
        ax.plot(xs, ys, "-x", label=utils.METHOD_LABEL[method],
                color=utils.METHOD_COLOR[method])
    if kwargs["title"]:
        ax.set_title("{energy_type} energy for {num_particles} particles "
                     "(ML = {ml}, K = {num_shells})"
                     .format(**locals()))
    ax.set_xlabel("ω")
    ax.set_ylabel("ε/ε[HF]")
    ax.legend()
    fig.tight_layout()
    utils.savefig(fig,
                  "by-freq-{num_shells}-{num_particles}-"
                  "{label}-{ml}-{interaction}"
                  .format(**locals()))
