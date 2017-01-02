import functools, json, os, sys
import matplotlib
import numpy as np
import pandas as pd
import scipy.optimize

JSON_PRETTY = {
    "ensure_ascii": False,
    "indent": 4,
    "separators": (",", ": "),
    "sort_keys": True,
}

GS_METHOD_COLOR = {
    "hf": "#841a0b",
    "mp2": "#39b237",
    "imsrg": "#1351c4",
}

METHOD_COLOR = {
    "qdpt": "#39b237",
    "eom": "#841a0b",
    "eom_quads": "#a825bc",
    "cc": "#1351c4",
}

def matplotlib_try_enable_deterministic_svgs():
    # we want deterministic SVGs, but this isn't supported until matplotlib 2.0
    try:
        matplotlib.rcParams["svg.hashsalt"] = ""
    except KeyError:
        sys.stderr.write("Warning: Your matplotlib is too old. "
                         "SVG output will be nondeterministic.\n")
        sys.stderr.flush()

def init(filename):
    os.chdir(os.path.dirname(filename))
    matplotlib_try_enable_deterministic_svgs()
    matplotlib.style.use("ggplot")

def skip_comment_char(read_func, filename):
    with open(filename) as f:
        s = f.read(2)
        assert s == "# "
        return read_func(f)

def load_json_records(fn):
    with open(fn) as f:
        return pd.DataFrame.from_records([
            json.loads(s) for s in f.read().split("\n\n") if s])

def sanitize_json(j):
    if isinstance(j, dict):
        return dict((k, sanitize_json(v)) for k, v in j.items())
    if isinstance(j, list) or isinstance(j, np.ndarray):
        return [sanitize_json(x) for x in j]
    if isinstance(j, float) or isinstance(j, int) or isinstance(j, str):
        return j
    if isinstance(j, np.float64):
        return float(j)
    if isinstance(j, np.int64):
        return int(j)
    raise TypeError("can't convert to JSON: {!r} ({})"
                    .format(j, type(j)))

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
    '''["num_shells", "num_filled", "freq", "ml", "label", "energy"]'''

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
        d["energy"] *= d["freq"] ** 0.5
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_IMSRG_FEI_HAM_particle_removed.dat",
                        delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom_f"
        d["energy"] *= d["freq"] ** 0.5
        d["interaction"] = v
        yield d

        # eom_quads: Nathan's EOM using Magnus method with quadruples

        d = pd.read_csv("EOM_magnus_quads_attached.dat",
                        delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "eom_quads"
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_magnus_quads_removed.dat",
                        delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom_quads"
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

def fit_change(fit_type, fit_points, x, y, **params):
    x_in = x[-fit_points:]
    y_in = y[-fit_points:]
    if fit_type == "loglog":
        def f(x, m, c):
            return np.exp(c) * x ** m
        p0 = np.polyfit(np.log(x_in), np.log(y_in), 1)
    else:
        def f(x, m, c):
            return np.exp(m * x + c)
        p0 = np.polyfit(x_in, np.log(y_in), 1)
    p, var_p = scipy.optimize.curve_fit(f, x_in, y_in, p0=p0)
    # chi-squared value per degree of freedom (using weight = 1)
    # normalized by
    badness = (np.sum((f(x_in, *p) - y_in) ** 2) / (len(x_in) - len(p)))
    r = {
        "fit_type": fit_type,
        "badness": float(badness),
        "params": p.tolist(),
        "cov_params": var_p,
    }
    r.update(params)
    if fit_type == "loglog":
        b = p[0] + 1.0
        a = np.exp(p[1]) / b
        jac_ba_mc = np.array([
            [1, 0],
            [-a / b, a],
        ])
        var_ba = jac_ba_mc.dot(var_p.dot(jac_ba_mc.transpose()))
        r["exponent"] = b
        r["exponent_err"] = var_ba[0, 0] ** 0.5
        r["coefficient"] = a
        r["coefficient_err"] = var_ba[1, 1] ** 0.5
        r["cov_exponent_coefficient"] = var_ba[0, 1]
    else:
        b = -1.0 / p[0]
        a = b * np.exp(p[1])
        jac_ba_mc = np.array([
            [b ** 2, 0],
            [b * a, a],
        ])
        var_ba = jac_ba_mc.dot(var_p.dot(jac_ba_mc.transpose()))
        r["lifetime"] = b
        r["lifetime_err"] = var_ba[0, 0] ** 0.5
        r["coefficient"] = a
        r["coefficient_err"] = var_ba[1, 1] ** 0.5
        r["cov_lifetime_coefficient"] = var_ba[0, 1]
    return {"x": x, "y": f(x, *p), "extra": r, "guess": p0}
