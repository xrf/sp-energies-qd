#!/usr/bin/env python3
import functools, os, sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import utils

utils.init(__file__)

d = utils.parse_simple("compare_rel_slopes.txt")

# we want to use (interaction, freq, num_particles) as an axis
# but matplotlib doesn't let us do that
# so we have to use this awful hack
x_labels = [
    "'_sigmaA=0.5_sigmaB=4.0', 1.0, 6",
    "'_sigmaA=0.5_sigmaB=4.0', 1.0, 12",
    "'_sigmaA=0.5_sigmaB=4.0', 1.0, 20",
    "'_sigmaA=0.5_sigmaB=4.0', 0.28, 6",
    "'_sigmaA=0.5_sigmaB=4.0', 0.28, 12",
    "'_sigmaA=0.5_sigmaB=4.0', 0.28, 20",
    "'_sigmaA=0.5', 1.0, 6",
    "'_sigmaA=0.5', 1.0, 12",
    "'_sigmaA=0.5', 1.0, 20",
    "'_sigmaA=0.5', 0.28, 6",
    "'_sigmaA=0.5', 0.28, 12",
    "'_sigmaA=0.5', 0.28, 20",
    "'_sigmaB=4.0', 1.0, 6",
    "'_sigmaB=4.0', 1.0, 12",
    "'_sigmaB=4.0', 1.0, 20",
    "'_sigmaB=4.0', 0.28, 6",
    "'_sigmaB=4.0', 0.28, 12",
    "'_sigmaB=4.0', 0.28, 20",
    "'', 1.0, 6",
    "'', 1.0, 12",
    "'', 1.0, 20",
    "'', 0.28, 6",
    "'', 0.28, 12",
    "'', 0.28, 20",
]
x_labels_ = [l.replace("_", " ").replace("sigma", "σ") for l in x_labels]
x_labels_lookup = dict((x, i) for i, x in enumerate(x_labels))
d["x_id"] = functools.reduce(lambda x, y: x + ", " + y, [
    d["interaction"],
    d["freq"].astype(str),
    d["num_particles"].astype(str),
]).map(lambda x: x_labels_lookup[x])

fig = plt.figure()
ax = fig.add_axes([0.1, 0.45, 0.8, 0.5])
for label, g in d.groupby(["label"]):
    ax.plot(g["x_id"],
            abs(g["slope"]),
            "o",
            markeredgewidth=0,
            label=label)
ax.set_ylabel("|rel_slope|")
ax.set_xlim(-1, len(x_labels))
ax.set_xticks(range(len(x_labels)))
ax.set_xticklabels(x_labels_, rotation=90)
ax.set_yscale("log")
ax.legend()
fn = "../FigureFiles/fig-rel-slopes.svg"
sys.stderr.write(fn + "\n")
sys.stderr.flush()
fig.savefig(fn)
plt.close(fig)
