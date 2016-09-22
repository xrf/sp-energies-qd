#!/usr/bin/env python3
import functools, os
import matplotlib.pyplot as plt
import pandas as pd

def skip_comment_char(read_func, filename):
    with open(filename) as f:
        s = f.read(2)
        assert s == "# "
        return read_func(f)

os.chdir(os.getcwd())

d = skip_comment_char(functools.partial(pd.read_csv, delim_whitespace=True),
                      "imsrg-qdpt/dat_arenergy_by_ml.txt")
d = d[["num_shells", "num_filled", "freq", "ml", "label", "energy"]]
d["method"] = "qdpt"
dq = d

d = pd.read_csv("EOM_IMSRG_qd_attached.dat", delim_whitespace=True)
d["energy"] = d["E(N+1)-E(N)"]
d = d[["shells", "filled", "ML", "omega", "energy"]]
d = d.rename(columns={
    "shells": "num_shells",
    "filled": "num_filled",
    "ML": "ml",
    "omega": "freq",
})
d["label"] = "add"
d["method"] = "eom"
dea = d

d = pd.read_csv("EOM_IMSRG_qd_removed.dat", delim_whitespace=True)
d["energy"] = -d["E(N-1)-E(N)"]
d = d[["shells", "filled", "ML", "omega", "energy"]]
d = d.rename(columns={
    "shells": "num_shells",
    "filled": "num_filled",
    "ML": "ml",
    "omega": "freq",
})
d["label"] = "rm"
d["method"] = "eom"
der = d

d = pd.concat([dq, dea, der])

num_filled = 2
freq = 0.28
label = "rm"

for num_filled in [2, 3]:
    for label in ["add", "rm"]:
        fig, ax = plt.subplots()
        ax.axhline(0.0, color="grey")
        for ml in [0, 1, 2]:
            dr = {
                "num_shells": [],
                "rel_diff": [],
            }
            for num_shells in [5, 6, 7, 8, 9, 10, 12, 13, 14, 15]:
                case = d[(d.num_shells == num_shells) &
                         (d.num_filled == num_filled) &
                         (d.freq == freq) &
                         (d.ml == ml) &
                         (d.label == label)]
                try:
                    energy_qdpt, = case[case.method == "qdpt"].energy
                    energy_eom, = case[case.method == "eom"].energy
                except ValueError:
                    continue
                rel_diff = ((energy_eom - energy_qdpt) /
                            (energy_eom + energy_qdpt) * 2.0)
                dr["num_shells"].append(num_shells)
                dr["rel_diff"].append(rel_diff)
            dr = pd.DataFrame.from_dict(dr)
            ax.plot(dr["num_shells"], dr["rel_diff"] * 100.0, "-x",
                    label="ml={}".format(ml))

        ax.set_title("number_of_particles = {}, type = {}"
                     .format(num_filled * (num_filled + 1),
                             {"add": "addition", "rm": "removal"}[label]))
        ax.set_xlabel("number_of_shells")
        ax.set_ylabel("relative_difference /%")
        ax.set_ylim(-1.5, 1.5)
        ax.legend()
plt.show()
