import collections, contextlib, functools, json, os, sys
import matplotlib
import matplotlib.pyplot as plt
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
    "eom_f": "#f4ad42",
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
    plot(filename).__enter__()

@contextlib.contextmanager
def plot(filename):
    os.chdir(os.path.dirname(filename))
    matplotlib_try_enable_deterministic_svgs()
    matplotlib.style.use("ggplot")
    yield
    if matplotlib.rcParams["interactive"]:
        plt.show(block=True)

def savefig(fig, name):
    if not matplotlib.rcParams["interactive"]:
        fn = "../FigureFiles/fig-{name}.svg".format(**locals())
        fig.savefig(fn)
        plt.close(fig)
        sys.stderr.write("// Figure saved to: {}\n\n".format(fn))
        sys.stderr.flush()

def groupby(d, cols):
    '''Like DataFrame.groupby but converts the key into an OrderedDict.'''
    for key, g in d.groupby(cols):
        # pandas inconsistently returns either a tuple or the object itself
        # depending on how many cols were used
        if len(cols) == 1:
            yield collections.OrderedDict([(cols[0], key)]), g
        else:
            yield collections.OrderedDict(zip(cols, key)), g

def update_range(r, x):
    for x in np.array(x).flatten():
        # make sure NaN gets handled correctly
        if r[0] == None or not (x >= r[0]):
            r[0] = x
        if r[1] == None or not (x <= r[1]):
            r[1] = x

def expand_range(r, factor):
    return np.array([
        r[0] - factor * (r[1] - r[0]),
        r[1] + factor * (r[1] - r[0]),
    ])

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

def json_pretty(j):
    return json.dumps(sanitize_json(j), **JSON_PRETTY)

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
    # note: y here is the derivative
    x_in = x[-fit_points:]
    y_in = y[-fit_points:]
    sign = np.sign(np.mean(y_in))
    if fit_type == "loglog":
        def f(x, b, a):
            return a * b * x ** (b - 1)
        mc0 = np.polyfit(np.log(x_in), np.log(abs(y_in)), 1)
        ba0 = [mc0[0] + 1.0, sign * np.exp(mc0[1]) / (mc0[0] + 1.0)]
        names = ["exponent", "coefficient"]
    elif fit_type == "semilog":
        def f(x, b, a):
            return -a * np.exp(-x / b) / b
        mc0 = np.polyfit(x_in, np.log(abs(y_in)), 1)
        ba0 = [-1.0 / mc0[0], sign * -np.exp(mc0[1]) / mc0[0]]
        names = ["lifetime", "coefficient"]
    else:
        assert False
    r0 = {
        "fit_type": fit_type + "0",
        "loglog_params": mc0,
    }
    for i, name in enumerate(names):
        r0[name] = ba0[i]

    # now try to do a nonlinear fit
    try:
        ba, var_ba = scipy.optimize.curve_fit(f, x_in, y_in, p0=ba0)
    except RuntimeError: # minimization failure
        return {"extra0": r0}
    # chi-squared value per degree of freedom (using weight = 1)
    badness = (np.sum((f(x_in, *ba) - y_in) ** 2) / (len(x_in) - len(ba)))
    r = {
        "fit_type": fit_type,
        "badness": badness,
    }
    r.update(params)
    for i, name in enumerate(names):
        r[name] = ba[i]
        r[name + "_err"] = var_ba[i, i] ** 0.5
        for j, name2 in enumerate(names[(i + 1):]):
            r["cov_{}_{}".format(name, name2)] = var_ba[i, j]

    return {"x": x, "y": f(x, *ba), "extra0": r0, "extra": r}
