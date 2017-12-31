import collections, functools, itertools, multiprocessing, sys
import numpy as np
import pandas as pd
import plot_fit, utils

DEFAULT_FIT_COUNT = 5
DEFAULT_MAXFEV = 100000
# fit_count=5:
#  maxfev=default: 198 fits failed, out of 2247
#  maxfev=10k: 40 fits failed, out of 2248
#  maxfev=100k: 0 fits failed

MIN_FIT_START = 4
MAX_FIT_START = 19

FIT_PARAM_KEYS = [
    "fit_stop",
    "coefficient",
    "coefficient_err",
    "exponent",
    "exponent_err",
    "constant",
    "constant_err",
]

def gather_fit_data(group, fit_count, maxfev):
    badness_threshold = plot_fit.LOGDERIV_BADNESS_THRESHOLD
    (label, interaction, num_filled, freq, method), gg = group
    context = collections.OrderedDict([
        ("label", label),
        ("interaction", interaction),
        ("num_filled", num_filled),
        ("freq", freq),
        ("method", method),
    ])
    gg = gg.rename(columns={"num_shells": "x", "energy": "y"})
    gg = gg.sort_values(["x"])
    deriv_gg = plot_fit.differentiate(gg, "x", "y", "dydx")
    results = []
    missing_num_shells = set()
    for fit_start in range(max(MIN_FIT_START, num_filled), MAX_FIT_START + 1):
        fit_stop = fit_start + fit_count - 1
        g = gg[gg["x"].between(fit_start, fit_stop)]
        if len(g) != fit_count:
            missing_num_shells.update(frozenset(range(fit_start, fit_stop + 1))
                                      - frozenset(gg["x"]))
            continue
        deriv_g = deriv_gg[deriv_gg["x"].between(fit_start, fit_stop)]
        fit = plot_fit.do_fit(g, deriv_g, badness_threshold=badness_threshold,
                              maxfev=maxfev)
        if fit is not None:
            main_fit = fit.get("full", fit["fixedab"])
            exponent = main_fit["exponent"]
            constant = main_fit["constant"]
            constant_err = main_fit["constant_err"]
            fit_info = {
                "fit_method": "full" if "full" in fit else "fixedab",
                "fit_stop": fit_stop,
                "chisq": np.linalg.norm(
                    g["y"]
                    - (main_fit["coefficient"] * g["x"] ** main_fit["exponent"]
                       + main_fit["constant"])) ** 2,
                "coefficient": main_fit["coefficient"],
                "coefficient_err": main_fit.get("coefficient_err", float("nan")),
                "exponent": exponent,
                "exponent_err": main_fit.get("exponent_err", float("nan")),
                "constant": constant,
                "constant_err": constant_err,
                "fixedab_constant_err": fit["fixedab"]["constant_err"],
                "dist": g["y"].tail(1).iloc[0] - constant,
                "dist_err": constant_err,
            }
        else:
            fit_info = {}
        results.append(fit_info)
    results = [utils.merge_dicts(context, x) for x in results]
    missing = utils.merge_dicts(context, {
        "num_shells": ",".join(map(str, sorted(missing_num_shells))),
    })
    return results, missing

def load_full_fit_data(fit_count=DEFAULT_FIT_COUNT,
                       maxfev=DEFAULT_MAXFEV):
    '''Load fit data from file if available.  Otherwise calculate the fits.'''
    fn = "fit_data.fit_count={fit_count}_maxfev={maxfev}.txt".format(**locals())
    try:
        return utils.load_table(fn)
    except OSError:
        pass

    sys.stderr.write("Fit data has not yet been calculated.  "
                     "This may take a few minutes...\n")
    sys.stderr.flush()
    d = utils.filter_preferred_ml(utils.load_all())
    d = d[~d["method"].isin(["imsrg[f]+eom[n]"])]
    with multiprocessing.Pool(4) as p:
        results_s, missing_num_shells = zip(*p.map(
            functools.partial(gather_fit_data,
                              fit_count=fit_count,
                              maxfev=maxfev),
            tuple(d.groupby(["label", "interaction", "num_filled",
                             "freq", "method"]))))
    results = itertools.chain(*results_s)

    missing_fn = ("fits_missing_points."
                  "fit_count={fit_count}_maxfev={maxfev}.log"
                  .format(**locals()))
    utils.save_table(missing_fn.format(**locals()),
                     pd.DataFrame.from_records(missing_num_shells))
    sys.stderr.write("Missing data points logged to: {}\n".format(missing_fn))
    sys.stderr.flush()

    d = pd.DataFrame.from_records(results)
    num_failed = (d["fit_method"] == "fixedab").sum()
    if num_failed:
        sys.stderr.write("{} out of {} fits failed\n"
                         .format(num_failed, len(d)))
        sys.stderr.flush()

    utils.save_table(fn, d)
    return d

def load_fit_data(label, fit_count=DEFAULT_FIT_COUNT, maxfev=DEFAULT_MAXFEV):
    d = load_full_fit_data(fit_count=fit_count, maxfev=maxfev)
    d = d[d["interaction"] == "normal"]
    d = d[d["label"] == label]
    d["num_particles"] = d["num_filled"] * (d["num_filled"] + 1)
    d = d.set_index(["fit_stop", "num_particles", "freq", "method"])
    return d

def get_fit_params(fit_data, num_shells, num_particles, freq, method):
    '''Raises KeyError if the data cannot be found.'''
    key = (num_shells, num_particles, freq, method)
    if key not in fit_data.index:
        raise KeyError()
    return utils.merge_dicts(dict(fit_data.loc[key]), {"fit_stop": num_shells})

def gather_predictive_data(group):
    (label, interaction, num_filled, freq, method), subresults = group
    subresults = [x for _, x in subresults.iterrows()]
    results = []
    for r in subresults[:-1]:
        best_constant = subresults[-1]["constant"]
        best_constant_err = subresults[-1]["constant_err"]
        constant = r["constant"]
        constant_err = r["constant_err"]
        dist = r["dist"]
        results.append(collections.OrderedDict([
            ("interaction", interaction),
            ("label", label),
            ("freq", freq),
            ("num_filled", num_filled),
            ("method", method),
            ("best_chisq", subresults[-1]["chisq"]),
            ("best_coefficient", subresults[-1]["coefficient"]),
            ("best_coefficient_err", subresults[-1]["coefficient_err"]),
            ("best_constant", best_constant),
            ("best_constant_err", best_constant_err),
            ("best_fixedab_constant_err", subresults[-1]["fixedab_constant_err"]),
            ("best_exponent", subresults[-1]["exponent"]),
            ("best_exponent_err", subresults[-1]["exponent_err"]),
            ("best_fit_method", subresults[-1]["fit_method"]),
            ("best_fit_stop", subresults[-1]["fit_stop"]),
            ("chisq", r["chisq"]),
            ("coefficient", r["coefficient"]),
            ("coefficient_err", r["coefficient_err"]),
            ("constant", constant),
            ("constant_err", constant_err),
            ("fixedab_constant_err", r["fixedab_constant_err"]),
            ("exponent", r["exponent"]),
            ("exponent_err", r["exponent_err"]),
            ("fit_method", r["fit_method"]),
            ("fit_stop", r["fit_stop"]),
            ("rel_discrep", constant / best_constant - 1.0),
            ("rel_discrep_err", (
                (constant_err / best_constant) ** 2 +
                (constant * best_constant_err / best_constant ** 2) ** 2
            ) ** 0.5),
            ("rel_dist", dist / best_constant),
            ("rel_dist_err", (
                (r["dist_err"] / best_constant) ** 2 +
                (dist * best_constant_err / best_constant ** 2) ** 2
            ) ** 0.5),
        ]))
    return results

def load_predictive_data(fit_count=DEFAULT_FIT_COUNT,
                         maxfev=DEFAULT_MAXFEV):
    d = load_full_fit_data(fit_count=fit_count, maxfev=maxfev)
    results = itertools.chain(*map(
        gather_predictive_data,
        tuple(d.groupby(["label", "interaction", "num_filled",
                         "freq", "method"]))))
    return pd.DataFrame.from_records(results)
