import collections, contextlib, decimal, functools, json, os, re, sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize

import time

# used for the 'float_precision' argument in pandas.read_csv
# I would rather use "round_trip", but it causes segfaults :/
PRECISION = "high"
PRECISION = "round_trip"

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
def plot(filename, block=True):
    os.chdir(os.path.dirname(filename))
    matplotlib_try_enable_deterministic_svgs()
    matplotlib.style.use("ggplot")
    yield
    if matplotlib.rcParams["interactive"]:
        plt.show(block=block)

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

def load_json_records(fn):
    with open(fn) as f:
        return pd.DataFrame.from_records([
            json.loads(s) for s in f.read().split("\n\n") if s])

def sanitize_json(j):
    if isinstance(j, dict):
        return dict((k, sanitize_json(v)) for k, v in j.items())
    if isinstance(j, list) or isinstance(j, tuple) or isinstance(j, np.ndarray):
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

def load_json(fn, fallback=None):
    if os.path.exists(fn):
        with open(fn) as f:
            return json.load(f)
    return fallback

def save_json(fn, j):
    with open(fn + ".tmp", "w") as f:
        f.write(json_pretty(j))
    os.rename(fn + ".tmp", fn)

def sync_axes_lims(ax, settings, save):
    setting = settings.get("xlim")
    if setting:
        ax.set_xlim(setting)
    setting = settings.get("ylim")
    if setting:
        ax.set_ylim(setting)
    def on_xlims_changed(ax):
        settings["xlim"] = ax.get_xlim()
        save()
    def on_ylims_changed(ax):
        settings["ylim"] = ax.get_ylim()
        save()
    ax.callbacks.connect("xlim_changed", on_xlims_changed)
    ax.callbacks.connect("ylim_changed", on_ylims_changed)

def parse_uncertainty(s):
    '''Parse a numeric string of the form "<number>(<uncertainty>)".'''
    value, uncertainty = re.match(r"([^(]+)\((\d+)\)$", s).groups()
    exponent = decimal.Decimal(value).as_tuple().exponent
    uncertainty =  decimal.Decimal((0, list(map(int, uncertainty)), exponent))
    return value, uncertainty

def parse_simple(fn):
    return pd.read_csv(fn, delim_whitespace=True, float_precision=PRECISION)

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

@np.vectorize
def get_num_filled(num_particles, max_shells=100):
    for k in range(max_shells):
        if num_particles == k * (k + 1):
            return k
    raise ValueError("could not find num_filled for num_particles={}"
                     .format(num_particles))

class FunDepViolationError(Exception):
    def __init__(self, in_cols, out_cols, violations, *args):
        self.in_cols = in_cols
        self.out_cols = out_cols
        self.violations = violations
        super().__init__(self, in_cols, out_cols, violations, *args)

    def __str__(self):
        return "Functional dependency does not hold: {!r} -> {!r}\n{!r}".format(
            self.in_cols, self.out_cols, self.violations)

def check_fun_dep(d, cols, tolers, combiner=lambda x: x):
    '''Check if a functional dependency constraint holds in the data.  That
    is, check if the output columns are unique (up to some tolerance) for each
    combination of input columns.  Then, a combined version is returned based
    on the `combiner` function, which by default is the identity function and
    therefore does nothing.  Raises `FunDepViolationError` when violations are
    found.

    For example, `check_fundep(d, ["x", "y"], {"z": 0.1, "w": 0.2}, combiner)`
    checks whether each `(z, w)` pair is uniquely determined by `(x, y)` with
    an allowable tolerance of `0.1` for `z` and `0.2` for `w`.

    `combiner` is a function `(DataFrame) -> DataFrame` that is fed to
    `GroupBy.apply` to combine duplicated data.
    '''
    tolers = collections.OrderedDict(tolers.items())
    cols = list(cols)
    out_cols = list(tolers.keys())
    gs = d.groupby(cols)

    # force NaN into errors using dropna(axis=1);
    # also, reset_index apparently modifies the GroupBy object
    std = d.groupby(cols).std(ddof=0).reset_index()
    mask = std.dropna(axis=1)[out_cols] > list(tolers.values())
    if mask.any().any():
        reduced_mask = functools.reduce(
            lambda x, y: x | y, map(lambda x: mask[x], out_cols))
        violations = pd.concat([gs.get_group(tuple(k))
                                for _, k in std[reduced_mask][cols].iterrows()])
        raise FunDepViolationError(cols, out_cols, violations)

    return gs.apply(combiner).reset_index(drop=True)

def load_gs_dmc_energies():
    '''["freq", "num_filled", "method", "energy"]'''
    d = pd.read_csv("gs-dmc-joergen.txt",
                    float_precision=PRECISION,
                    header=0, index_col=False,
                    delim_whitespace=True, comment="#")
    d["num_filled"] = get_num_filled(d["num_particles"])
    del d["num_particles"]
    return d

def priority_combiner(d):
    d = d.drop_duplicates()
    d = d[d["priority"] == d["priority"].max()]
    if len(d) > 1:
        raise Exception("cannot resolve data with equal priority")
    return d

SAM_ATTACHED_COLS = ["shells", "filled", "ML", "MS", "omega", "E(N)",
                     "E(N+1)-E(N)", "E(N+1)", "partialnorm(1p)"]

SAM_REMOVED_COLS = ["shells", "filled", "ML", "MS", "omega", "E(N)",
                    "E(N-1)-E(N)", "E(N+1)", "partialnorm(1p)"]

def load_gs_energies():
    '''["freq", "num_filled", "num_shells", "method", "energy"]'''
    ds = []

    # agreement between Fei's IMSRG and Sarah's IMSRG ~ 1e-3
    # agreement between Nathan's IMSRG and Sarah's IMSRG ~ 1e-3
    d = pd.read_csv("../compWithOtherMethods/data_formatted.txt",
                    float_precision=PRECISION,
                    index_col=False, delim_whitespace=True, comment="#")
    # this data point does not agree with Fei's results and just seems wrong
    # (it breaks the monotonically decreasing behavior of HF)
    d = d[~((d["freq"] == 0.1) &
            (d["num_particles"] == 30) &
            (d["num_shells"] == 16) &
            (d["method"] == "hf"))]
    # this data point disagrees with both Fei and Nathan's results
    d = d[~((d["freq"] == 0.5) &
            (d["num_particles"] == 30) &
            (d["num_shells"] == 11) &
            (d["method"] == "imsrg"))]
    d["num_filled"] = get_num_filled(d["num_particles"])
    d["priority"] = -2
    del d["num_particles"]
    d = d.dropna()
    ds.append(d)

    d = pd.read_csv("imsrg-qdpt/dat_gsenergy.txt",
                    float_precision=PRECISION,
                    header=0, index_col=False,
                    delim_whitespace=True, comment="#")
    d["priority"] = 0
    del d["time"]
    ds.append(d)

    # agreement between Fei's IMSRG and Nathan's IMSRG ~ 1e-4
    d1 = pd.read_csv("EOM_IMSRG_qd_attached.dat",
                     float_precision=PRECISION,
                     delim_whitespace=True)
    d1["method"] = "imsrg"
    d2 = pd.read_csv("EOM_CCSD_qd_attached.dat",
                     float_precision=PRECISION,
                     header=None, names=SAM_ATTACHED_COLS,
                     delim_whitespace=True)
    d2["method"] = "cc"
    d = pd.concat([d1, d2], ignore_index=True)
    d = d[["shells", "filled", "omega", "E(N)", "method"]]
    d = d.rename(columns={
        "shells": "num_shells",
        "filled": "num_filled",
        "omega": "freq",
        "E(N)": "energy",
    })
    d["priority"] = -1
    d = d.drop_duplicates()
    ds.append(d)

    d = pd.concat(ds, ignore_index=True)
    try:
        d = check_fun_dep(d, ["freq", "num_filled", "num_shells", "method"],
                          {"energy": 6e-4}, combiner=priority_combiner)
    except FunDepViolationError as e:
        with open("fun_dep_violations.out", "w") as f:
            e.violations.to_csv(f, sep=" ", index=False)
        raise e

    del d["priority"]
    return d

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

        d = pd.read_csv("EOM_IMSRG_qd_attached.dat",
                        float_precision=PRECISION, delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "eom"
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_IMSRG_qd_removed.dat",
                        float_precision=PRECISION, delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom"
        d["interaction"] = v
        yield d

        d = pd.read_csv("freq_sweep_N6_R10_attached.dat",
                        float_precision=PRECISION, delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "eom"
        d["interaction"] = v
        yield d

        d = pd.read_csv("freq_sweep_N6_R10_removed.dat",
                        float_precision=PRECISION, delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom"
        d["interaction"] = v
        yield d

        # eom_f: Nathan's EOM on Fei's IMSRG matrix elements

        d = pd.read_csv("EOM_IMSRG_FEI_HAM_particle_attached.dat",
                        float_precision=PRECISION, delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "eom_f"
        d["energy"] *= d["freq"] ** 0.5
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_IMSRG_FEI_HAM_particle_removed.dat",
                        float_precision=PRECISION, delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom_f"
        d["energy"] *= d["freq"] ** 0.5
        d["interaction"] = v
        yield d

        # eom_quads: Nathan's EOM using Magnus method with quadruples

        d = pd.read_csv("EOM_magnus_quads_attached.dat",
                        float_precision=PRECISION, delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "eom_quads"
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_magnus_quads_removed.dat",
                        float_precision=PRECISION, delim_whitespace=True)
        d = parse_nathan_like_data(d, "rm")
        d["method"] = "eom_quads"
        d["interaction"] = v
        yield d

        # cc: Sam's coupled cluster

        d = pd.read_csv("EOM_CCSD_qd_attached.dat",
                        float_precision=PRECISION,
                        header=None,
                        names=NATHAN_ATTACHED_COLS,
                        delim_whitespace=True)
        d = parse_nathan_like_data(d, "add")
        d["method"] = "cc"
        d["interaction"] = v
        yield d

        d = pd.read_csv("EOM_CCSD_qd_removed.dat",
                        float_precision=PRECISION,
                        header=None,
                        names=NATHAN_REMOVED_COLS,
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
