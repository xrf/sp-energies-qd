#!/usr/bin/env python3
import argparse, itertools, os, sys
import matplotlib
import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import fits, utils

def plot(label, num_filled, freq,
         num_shells_range):
    fit_data = fits.load_fit_data(label)

    d = utils.load_all()
    d = utils.filter_preferred_ml(d)
    d = d[~d["method"].isin(["imsrg[f]+eom[n]",
                             "imsrg",
                             "magnus_quads+eom",
                             "fci",
                             "hf"]) &
          (d["interaction"] == "normal") &
          (d["label"] == label) &
          (d["num_filled"] == num_filled) &
          (d["num_shells"] >= num_shells_range[0]) &
          (d["num_shells"] <= num_shells_range[1]) &
          (d["freq"] == freq)]
    num_particles = num_filled * (num_filled + 1)
    energy_type = {"ground": "ground state",
                   "add": "addition",
                   "rm": "removal"}[label]
    fig, ax = plt.subplots()
    fig.set_size_inches((3.2, 2.7))
    base_markersize = 5
    xc = np.linspace(num_shells_range[0] - 0.5, num_shells_range[1] + 0.5, 200)
    for method, case in d.groupby("method"):
        case = case.sort_values("num_shells")
        xs = case["num_shells"].astype(int)
        ys = case["energy"]
        marker = utils.METHOD_MARKER[method]
        style = {
            "marker": marker,
            "markerfacecolor": "none",
            "markersize": (utils.MARKERSIZE_CORRECTION.get(marker, 1.0) *
                           base_markersize),
            "color": utils.METHOD_COLOR[method],
            "label": utils.METHOD_LABEL[method],
        }
        ax.plot(xs, ys, linestyle="none", **style)
        p = fits.get_fit_params(fit_data, np.max(xs), num_particles,
                                freq, method)
        ax.plot(xc, p["coefficient"] * xc ** p["exponent"] + p["constant"],
                **utils.merge_dicts(style, {"marker": ""}))
    ax.axvspan(np.min(xc), np.max(xs) - fits.DEFAULT_FIT_COUNT + 0.5,
               color="#cccccc")
    ax.set_xlim(np.min(xc), np.max(xc))
    ax.get_xaxis().set_major_locator(
        matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_xlabel("$K$ (number of shells)")
    ax.set_ylabel("$E$ (energy)")
    ax.legend(frameon=False)
    fig.tight_layout()
    utils.savefig(fig,
                  "../Manuscript/fig-fit-{num_filled}-{freq}-{label}.pdf"
                  .format(**locals()))

with utils.plot(__file__, call=plot) as interactive:
    if not interactive:
        plot("ground", freq=1.0, num_filled=2, num_shells_range=[6, 15])
