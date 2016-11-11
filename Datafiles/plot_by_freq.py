#!/usr/bin/env python3
# Generates ../FigureFiles/fig-by-freq-*.svg from the data files.
import itertools, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import utils

os.chdir(os.path.dirname(__file__))

d = pd.concat(utils.get_ar_energies())
print("Saving ...".format(**locals()), flush=True)
num_shells = 10.0
num_filled = 2
d = d[(d["num_shells"] == num_shells) &
      (d["num_filled"] == num_filled)]
num_particles = num_filled * (num_filled + 1)
for interaction, label in itertools.product(
        [utils.V0, utils.V2],
        ["add", "rm"]):
    energy_type = {"add": "addition", "rm": "removal"}[label]
    ml = (num_filled + (label != "add")) % 2
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    for method in d["method"].unique():
        case = d[(d["ml"] == ml) &
                 (d["label"] == label) &
                 (d["method"] == method) &
                 (d["interaction"] == interaction)]
        case = case.sort_values("freq").drop_duplicates()
        xs = case["freq"]
        ys = case["energy"]
        kwargs = {
            "qdpt": {
                "color": "#39b237",
            },
            "eom": {
                "color": "#841a0b",
            },
            "eom_f": {
                "color": "#c6cc26",
                "linestyle": "--",
            },
            "cc": {
                "color": "#1351c4",
            },
        }[method]
        ax.plot(xs, ys, "-x", label=method,
                linewidth=1.5, **kwargs)
    ax.set_title("{energy_type} energy for {num_particles} particles "
                 "(ML = {ml}, num_shells = {num_shells})"
                 .format(**locals()))
    ax.set_xlabel("freq")
    ax.set_ylabel("energy")
    ax.legend()
    fn = ("../FigureFiles/fig-by-freq-"
          "{num_shells}-{num_particles}-{label}-{ml}{interaction}.svg"
          .format(**locals()))
    sys.stderr.write(fn + "\n")
    sys.stderr.flush()
    fig.savefig(fn)
    plt.close(fig)
