import functools
from numpy import sqrt
import pandas as pd

def skip_comment_char(read_func, filename):
    with open(filename) as f:
        s = f.read(2)
        assert s == "# "
        return read_func(f)

def parse_simple(fn):
    return skip_comment_char(
        functools.partial(pd.read_csv, delim_whitespace=True), fn)

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

# interactions
V0 = ""                                 # normal Coulomb
V2 = "_sigmaA=0.5_sigmaB=4.0"           # softened Coulomb (see figures.md)

def get_ar_energies():
    yield from get_ar_energies_for_v(V0)
    yield from get_ar_energies_for_v(V2)

def get_ar_energies_for_v(v):

    # qdpt: Fei's QDPT on Fei's IMSRG matrix elements

    d = parse_simple("imsrg-qdpt/dat_arenergy_by_ml{v}.txt".format(**locals()))
    d = d[["num_shells", "num_filled", "freq", "ml", "label", "energy"]]
    d["method"] = "qdpt"
    d["interaction"] = v
    yield d

    if v == V0:

        # eom: Nathan's EOM on Nathan's IMSRG matrix elements

        d = pd.read_csv("EOM_IMSRG_qd_attached.dat", delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "eom"
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_IMSRG_qd_removed.dat", delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom"
        d["interaction"] = v
        yield d

        d = pd.read_csv("freq_sweep_N6_R10_attached.dat", delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "eom"
        d["interaction"] = v
        yield d

        d = pd.read_csv("freq_sweep_N6_R10_removed.dat", delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom"
        d["interaction"] = v
        yield d

        # eom_f: Nathan's EOM on Fei's IMSRG matrix elements

        d = pd.read_csv("EOM_IMSRG_FEI_HAM_particle_attached.dat",
                        delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "eom_f"
        d["energy"] *= sqrt(d["freq"])
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_IMSRG_FEI_HAM_particle_removed.dat",
                        delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom_f"
        d["energy"] *= sqrt(d["freq"])
        d["interaction"] = v
        yield d

        # cc: Sam's coupled cluster

        d = pd.read_csv("EOM_CCSD_qd_attached.dat", header=None,
                        names=["shells", "filled", "ML", "MS", "omega", "E(N)",
                               "E(N+1)-E(N)", "E(N+1)", "partialnorm(1p)"],
                        delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "cc"
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_CCSD_qd_removed.dat", header=None,
                        names=["shells", "filled", "ML", "MS", "omega", "E(N)",
                               "E(N-1)-E(N)", "E(N+1)", "partialnorm(1p)"],
                        delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "cc"
        d["interaction"] = v
        yield d
