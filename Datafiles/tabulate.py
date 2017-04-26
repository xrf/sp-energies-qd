#!/usr/bin/env python3
import itertools, math, os, sys
import numpy as np
import pandas as pd
import fits, utils

NUM_SHELLS = [12, 15, 15, 15, 20, 20]

# known cases of IM-SRG that don't converge
KNOWN_NC_IMSRG = set((k, kf * (kf + 1), f) for (k, kf, f) in [
    (12, 5, 0.1),
    (12, 6, 0.1),
    (12, 6, 0.1),
    (13, 6, 0.1),
    (15, 6, 0.1),
    (16, 7, 0.28),
    (20, 7, 0.1),
])

def render_entry(value, err=None, default_precision=4):
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
                            fit_count=fits.DEFAULT_FIT_COUNT,
                            maxfev=fits.DEFAULT_MAXFEV):
    d = fits.load_fit_data(label, fit_count=fit_count, maxfev=maxfev)

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
                        p = fits.get_fit_params(d, num_shells, num_particles,
                                                freq, method)
                    except fits.NoMatchButBestExistsError as e:
                        if "--verbose" in sys.argv:
                            sys.stderr.write("warning: {!r}\n".format(e))
                            sys.stderr.flush()
                        return np.nan,
                    except KeyError:
                        return np.nan,

                    value = p["constant"]
                    err = p["constant_err"]
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
