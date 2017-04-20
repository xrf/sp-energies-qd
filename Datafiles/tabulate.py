#!/usr/bin/env python3
import itertools, math, os, sys
import numpy as np
import pandas as pd
import utils

NUM_SHELLS = [12, 15, 15, 15, 20, 20]

# known cases of IM-SRG that don't converge
KNOWN_NC_IMSRG = set((k, kf * (kf + 1), f) for (k, kf, f) in [
    (20, 7, 0.1),
])

def render_entry(value, err=None, default_precision=5):
    if isinstance(value, str):
        return "{" + value + "}"
    if np.isnan(value):
        return ""
    if err is None:
        return "{{:.{}f}}".format(default_precision).format(value)
    if np.isnan(err):
        return "{{:.{}f}}{{{{(?)}}}}".format(default_precision).format(value)
    return utils.render_with_err(value, err)

def save_table(path, label,
               num_shells=NUM_SHELLS,
               num_particles=[6, 12, 20, 30, 42, 56],
               freqs=[0.1, 0.28, 1.0]):
    d = utils.load_all()
    d = utils.filter_preferred_ml(d)
    del d["ml"]
    d = d[d["interaction"] == "normal"]
    del d["interaction"]
    d = d[d["label"] == label]
    del d["label"]
    del d["num_filled"]
    d = d.set_index(["num_shells", "num_particles", "freq", "method"])

    if label == "ground":
        methods = ["hf", "mp2", "imsrg", "ccsd"]
        header = r"""
        \begin{tabular}{S[table-format=2.0]SS[table-format=2.0]S[table-format=3.5]S[table-format=3.5]S[table-format=3.5]S[table-format=3.5]}%
        \toprule
        {$N$} & {$\omega$} & {$K$} & {HF} & {MP2} & {IM-SRG(2)} & {CCSD} \\
        """
    else:
        methods = ["hf+qdpt3", "imsrg+qdpt3", "imsrg+eom", "ccsd+eom"]
        header = r"""
        \begin{tabular}{S[table-format=2.0]SS[table-format=2.0]S[table-format=3.5]S[table-format=3.5]S[table-format=3.5]S[table-format=3.5]}%
        \toprule
        {$N$} & {$\omega$} & {$K$} & {HF} & {IM-SRG(2)} & {IMSRG(2)} & {CCSD} \\
        {} & {} & {} & {+QDPT3} & {+QDPT3} & {+EOM} & {+EOM} \\
        """
    c = d["energy"]

    s = []
    s.append(header)
    for num_shells, num_particles in zip(num_shells, num_particles):
        s.append("\\midrule\n")
        for freq in freqs:
            row = [num_particles, freq, num_shells]
            for method in methods:
                def get(): # to allow early returns
                    if ("imsrg" in method and
                            (num_shells, num_particles, freq)
                                in KNOWN_NC_IMSRG):
                        return "n.c."
                    try:
                        r = d.loc[(num_shells, num_particles, freq, method)]
                    except KeyError:
                        return np.nan
                    return r["energy"]
                row.append(render_entry(get()))
            s.append(" & ".join(map(str, row)) + " \\\\\n")
    s.append(r"""\bottomrule\end{tabular}""")

    with open(path, "w") as f:
        f.write("".join(s))
    return d

def dump_best_fit_stop(d):
    d = d[d["freq"].isin([0.1, 0.28, 1.0])]
    d = d.groupby(["num_particles", "freq", "method"]).first()
    d = d[["best_fit_stop"]]
    d = d.unstack("method")
    sys.stdout.write(d)
    sys.exit(1)

def save_extrapolated_table(path, label,
                            num_shells=NUM_SHELLS,
                            num_particles=[6, 12, 20, 30, 42, 56],
                            freqs=[0.1, 0.28, 1.0],
                            fit_count=5,
                            maxfev=100000):
    d = utils.load_table(f"fits.fit_count={fit_count}_maxfev={maxfev}.txt")
    d = d[d["interaction"] == "normal"]
    d = d[d["label"] == label]
    d["num_particles"] = d["num_filled"] * (d["num_filled"] + 1)
    #dump_best_fit_stop(d)
    d2 = d.groupby(["num_particles", "freq", "method"]).first()
    d = d.set_index(["fit_stop", "num_particles", "freq", "method"])

    # HF has been excluded because:
    # (1) the extrapolations are kinda meh (because HF isn't very power law)
    # (2) the tables are already too wide
    if label == "ground":
        methods = ["mp2", "imsrg", "ccsd"]
        header = r"""
        \begin{tabular}{S[table-format=2.0]SS[table-format=2.0]S[table-format=4.6]S[table-format=4.6]S[table-format=4.6]}%
        \toprule
        {$N$} & {$\omega$} & {$K_{\text{stop}}$} & {MP2} & {IM-SRG(2)} & {CCSD} \\
        """
    else:
        methods = ["imsrg+qdpt3", "imsrg+eom", "ccsd+eom"]
        header = r"""
        \begin{tabular}{S[table-format=2.0]SS[table-format=2.0]S[table-format=3.6]S[table-format=3.6]S[table-format=3.6]}%
        \toprule
        {$N$} & {$\omega$} & {$K_{\text{stop}}$} & {IM-SRG(2)} & {IMSRG(2)} & {CCSD} \\
        {} & {} & {} & {+QDPT3} & {+EOM} & {+EOM} \\
        """

    s = []
    s.append(header)
    for num_shells, num_particles in zip(num_shells, num_particles):
        s.append("\\midrule\n")
        for freq in freqs:
            row = [num_particles, freq, num_shells]
            for method in methods:
                def get(): # to allow early returns

                    # check for nonconvergent cases
                    if method in ["imsrg+qdpt3"]:
                        for k in range(num_shells + 1 - fit_count,
                                       num_shells + 1):
                            if (k, num_particles, freq) in KNOWN_NC_IMSRG:
                                return "n.c.",

                    try:
                        r = d.loc[(num_shells, num_particles, freq, method)]
                    except KeyError:
                        r = None
                    else:
                        value = r["constant"]
                        err = r["constant_err"]

                    if r is None:
                        # some of the fit_stops are "best"
                        # so they aren't included directly
                        try:
                            r = d2.loc[(num_particles, freq, method)]
                        except KeyError:
                            return np.nan,
                        if r["best_fit_stop"] != num_shells:
                            sys.stderr.write(
                                "no matching fit "
                                "even though best exists: {} -> {}\n"
                                .format((num_particles, freq, method),
                                        r["best_fit_stop"]))
                            return np.nan,
                        value = r["best_constant"]
                        err = r["best_constant_err"]

                    if err > abs(value):
                        value = "{n.f.}"
                    return value, err

                row.append(render_entry(*get()))
            s.append(" & ".join(map(str, row)) + " \\\\\n")
    s.append(r"""\bottomrule\end{tabular}""")

    with open(path, "w") as f:
        f.write("".join(s))
    return d

os.chdir(os.path.dirname(__file__))
save_table("../Manuscript/tab-ground.tex", label="ground")
save_table("../Manuscript/tab-add.tex", label="add")
save_table("../Manuscript/tab-rm.tex", label="rm")
save_extrapolated_table("../Manuscript/tab-ground-extrapolated.tex",
                        label="ground")
save_extrapolated_table("../Manuscript/tab-add-extrapolated.tex", label="add")
save_extrapolated_table("../Manuscript/tab-rm-extrapolated.tex", label="rm")
