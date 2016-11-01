#!/usr/bin/env python3
# Generates ../FigureFiles/fig-compare-*.svg from the data files.
import functools, os
import matplotlib.pyplot as plt
import pandas as pd
from numpy import sqrt
import utils

def parse_nathan_like_data(d, label):
    if label == "add":
        d["energy"] = d["E(N+1)-E(N)"]
    elif label == "rm":
        d["energy"] = -d["E(N-1)-E(N)"]
    else:
        raise ValueError("invalid value for label parameter")
    d = d[["shells", "filled", "ML", "omega", "energy"]]
    d = d.rename(columns={
        "shells": "num_shells",
        "filled": "num_filled",
        "ML": "ml",
        "omega": "freq",
    })
    d["method"] = "eom"
    d["label"] = label
    return d

os.chdir(os.path.dirname(__file__))

ds = []

d = utils.skip_comment_char(
    functools.partial(pd.read_csv, delim_whitespace=True),
    "imsrg-qdpt/dat_arenergy_by_ml.txt")
d = d[["num_shells", "num_filled", "freq", "ml", "label", "energy"]]
d["method"] = "qdpt"
ds.append(d)

d = pd.read_csv("EOM_IMSRG_qd_attached.dat", delim_whitespace=True)
d = parse_nathan_like_data(d, "add")
d["method"] = "eom"
ds.append(d)

d = pd.read_csv("EOM_IMSRG_qd_removed.dat", delim_whitespace=True)
d = parse_nathan_like_data(d, "rm")
d["method"] = "eom"
ds.append(d)

d = pd.read_csv("EOM_IMSRG_FEI_HAM_particle_attached.dat", delim_whitespace=True)
d = parse_nathan_like_data(d, "add")
d["method"] = "eom_f"
d["energy"] *= sqrt(d["freq"])
ds.append(d)

d = pd.read_csv("EOM_IMSRG_FEI_HAM_particle_removed.dat", delim_whitespace=True)
d = parse_nathan_like_data(d, "rm")
d["method"] = "eom_f"
d["energy"] *= sqrt(d["freq"])
ds.append(d)

d = pd.read_csv("EOM_CCSD_qd_attached.dat", header=None,
                names=["shells", "filled", "ML", "MS", "omega",
                       "E(N)", "E(N+1)-E(N)", "E(N+1)", "partialnorm(1p)"],
                delim_whitespace=True)
d = parse_nathan_like_data(d, "add")
d["method"] = "cc"
ds.append(d)

d = pd.read_csv("EOM_CCSD_qd_removed.dat", header=None,
                names=["shells", "filled", "ML", "MS", "omega",
                       "E(N)", "E(N-1)-E(N)", "E(N+1)", "partialnorm(1p)"],
                delim_whitespace=True)
d = parse_nathan_like_data(d, "rm")
d["method"] = "cc"
ds.append(d)

d = pd.concat(ds)

print("Saving ...".format(**locals()), flush=True)
for freq in [0.28, 1.0]:
    for num_filled in [2, 3, 4]:
        num_particles = num_filled * (num_filled + 1)
        for label in ["add", "rm"]:
            fig, axs = plt.subplots(2, 1)
            fig.set_size_inches(8, 10)
            ml = (num_filled + (2 != "add")) % 2
            for zoomed in [False, True]:
                ax = axs[int(zoomed)]
                # zoom in to the part where it's nearly converged
                num_shells_start = {
                    2: 7,
                    3: 9,
                    4: 11,
                }[num_filled] if zoomed else 5
                for method in d["method"].unique():
                    xs = []
                    ys = []
                    for num_shells in range(num_shells_start, 16):
                        case = d[(d["num_shells"] == num_shells) &
                                 (d["num_filled"] == num_filled) &
                                 (d["freq"] == freq) &
                                 (d["ml"] == ml) &
                                 (d["label"] == label)]
                        try:
                            energy, = case[case["method"] == method].energy
                        except ValueError:
                            continue
                        xs.append(num_shells)
                        ys.append(energy)
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
                if not zoomed:
                    ax.set_title("{} energy for {} particles (ML = {}, omega = {})"
                                 .format({"add": "addition",
                                          "rm": "removal"}[label],
                                         num_particles, ml, freq))
                    ax.legend()
                ax.set_xlabel("number_of_shells")
                ax.set_ylabel("energy")
            fn = ("../FigureFiles/fig-compare-"
                  "{num_particles}-{freq}-{label}-{ml}.svg"
                  .format(**locals()))
            print(fn, flush=True)
            fig.savefig(fn)
            plt.close(fig)
