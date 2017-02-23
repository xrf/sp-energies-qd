#!/usr/bin/env python3
import itertools, os, sys
import numpy as np
import pandas as pd
import utils

def save_table(path, label,
               num_shells=[15, 15, 15, 15, 15, 20],
               num_particles=[6, 12, 20, 30, 42, 56],
               freqs=[0.1, 0.28, 1.0]):
    d = utils.load_all()
    d = utils.filter_preferred_ml(d)
    print(d["method"].unique())
    del d["ml"]
    d = d[d["interaction"] == "normal"]
    del d["interaction"]
    d = d[d["label"] == label]
    del d["label"]
    del d["num_filled"]
    d = d.set_index(["num_shells", "num_particles", "freq", "method"])
    d = d.unstack("method")

    if label == "ground":
        methods = ["hf", "mp2", "imsrg", "ccsd"]
        header = r"""
        \begin{tabular}{SSSSSSS}%
        \toprule
        {$N$} & {$\omega$} & {$K$} & {HF} & {MP2} & {IM-SRG(2)} & {CCSD} \\
        """
    else:
        methods = ["hf+qdpt3", "imsrg+qdpt3", "imsrg+eom", "ccsd+eom"]
        header = r"""
        \begin{tabular}{SSSSSSS}%
        \toprule
        {$N$} & {$\omega$} & {$K$} & {HF} & {IM-SRG(2)} & {IMSRG(2)} & {CCSD} \\
        {} & {} & {} & {+QDPT3} & {+QDPT3} & {+EOM} & {+EOM} \\
        """
    d = d.reindex_axis(methods, level=-1, axis=1)
    d.columns = d.columns.droplevel(level=0)

    out_file = open(path, "w")
    out_file.write(header)
    for num_shells, num_particles in zip(num_shells, num_particles):
        out_file.write("\\midrule\n")
        for freq in freqs:
            row = []
            row.extend([num_particles, freq, num_shells])
            for method in methods:
                value = d.loc[(num_shells, num_particles, freq), method]
                row.append("" if np.isnan(value) else
                           "{:.4f}".format(value))
            out_file.write(" & ".join(map(str, row)) + " \\\\\n")
    out_file.write(r"""\bottomrule\end{tabular}""")

    return d

os.chdir(os.path.dirname(__file__))
save_table("../Manuscript/tab-ground.tex", label="ground")
save_table("../Manuscript/tab-add.tex", label="add")
save_table("../Manuscript/tab-rm.tex", label="rm")
