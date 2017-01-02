#!/usr/bin/env python3
# Generates ../FigureFiles/fig-fit-badness.svg from the data files.
import argparse, json, os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utils

FIT_TYPES = ["loglog", "semilog"]

def get_group_or(d, group, key, default=float("nan")):
    try:
        g = d.get_group(group)
    except KeyError:
        return default
    if len(g) == 0:
        return default
    x, = g[key]
    return x

def markersize(marker):
    return {
        "*": 1.2,
    }.get(marker, 1.0)

def plot(interactive):
    if interactive:
        plt.ion()

    d = utils.load_json_records("fits.txt")

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.45, 0.6, 0.5])
    fig.set_size_inches(8, 10)

    x_ticks = sorted(d["title"].unique())
    x_indices = range(len(x_ticks))
    for [fit_type, method], g in d.groupby(["fit_type", "method"]):
        s = g.groupby("title")
        marker = {
            "qdpt": "o",
            "eom": "x",
            "eom_quads": "v",
            "cc": "*",
        }[method]
        ax.plot(
            x_indices,
            [get_group_or(s, x, "badness") for x in x_ticks],
            "o",
            label="{}, {}".format(fit_type, method),
            markersize=(10.0 * markersize(marker)),
            marker=marker,
            color={
                "semilog": "#d11987",
                "loglog": "#088fdd",
            }[fit_type],
        )
    ax.legend(bbox_to_anchor=(1.45, 0.8))
    ax.set_xticks(x_indices)
    ax.set_xticklabels(x_ticks, rotation=90)
    ax.set_ylabel("badness")
    ax.set_yscale("log")

    if not interactive:
        fn = "fig-fit-badness.svg"
        fig.savefig(fn)
        plt.close(fig)
        sys.stderr.write("// Figure saved to: {}\n\n".format(fn))
        sys.stderr.flush()

    if interactive:
        plt.show(block=True)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--interactive", action="store_true")
    kwargs = vars(p.parse_args())
    utils.init(__file__)
    plot(**kwargs)

main()
