#!/usr/bin/env python3
import argparse, itertools, json, logging, os, sys, traceback, warnings
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
import scipy, scipy.stats
import utils

DUPLICATE_THRESHOLD = 1e-6

LOGDERIV_BADNESS_THRESHOLD = 0.03

STAGE_TO_LINESTYLE = {
    "logderiv": ":",
    "fixedab": "--",
    "full": "-",
}

def default_fit_condition(cov, **kwargs):
    return np.isfinite(cov).all()

def rel_badness_cov(ref_fit, fit):
    # dimensionally valid and has the advantage of not being dependent
    # on the parameter values, but I'm not sure if there is a rigorous
    # justification to its sensibility
    n = ref_fit["cov"].shape[0]
    with np.errstate(divide="ignore"):
        return np.linalg.norm(fit["cov"] / ref_fit["cov"]) / n ** 2

def try_curve_fit(f, x, y, p0, badness_threshold, min_num_points=None,
                  rel_badness_func=rel_badness_cov, max_num_outliers=2):
    '''Like scipy.optimize.curve_fit, but tries to drop some of the points to
    find the most satisfactory fit.

    rel_badness_func(ref_fit, fit) should return a number between 0.0 and inf,
    where 0.0 = fit is better than ref_fit, 1.0 = equally good, etc.  A 'fit'
    is a dictionary:

        {"outliers": outlier_indices, "p": params, "cov": covariance}

    This is also the value that gets returned.

    The badness_threshold is a number between 0.0 and 1.0 that determines the
    eagerness to drop points.  At 1.0 the algorithm is maximally eager.
    '''
    assert 0.0 <= badness_threshold <= 1.0
    if min_num_points is None:
        min_num_points = len(p0) + 2
    max_num_outliers = max(min(max_num_outliers, len(x) - min_num_points), 0)
    ref_fit = None
    best_fit = None
    threshold = badness_threshold
    for num_outliers in range(max_num_outliers + 1):
        fits = []
        for outliers in itertools.combinations(list(x.index), num_outliers):
            outliers = list(outliers)
            new_x = x.drop(outliers)
            new_y = y.drop(outliers)
            try:
                p, cov = scipy.optimize.curve_fit(f, new_x, new_y, p0=p0)
            except TypeError:
                return
            fit = {"outliers": outliers, "p": p, "cov": cov}
            if ref_fit is None:         # first fit (no excluded points)
                ref_fit = fit
                best_fit = fit
            else:
                fits.append(fit)
        if fits:
            rel_badnesses = [rel_badness_func(ref_fit, fit) for fit in fits]
            change = (np.min(rel_badnesses)
                      / scipy.stats.gmean(rel_badnesses)
                      / threshold)
            if change <= 1.0:
                threshold *= badness_threshold
                best_fit = fits[np.argmin(rel_badnesses)]
                logging.info("try_curve_fit: accept: "
                             "log10(min_badness/avg_badness) = "
                             f"{np.log10(change * badness_threshold)}")
        elif badness_threshold == 0.0:
            break
    return best_fit

def do_fit(data, deriv_data, badness_threshold, maxfev=0):
    '''Perform a fit using a power law model:

        y = a * x ** b + c

    data: DataFrameish{"x", "y"}
    deriv_data: DataFrameish{"x", "dydx"}
    '''

    result = {}

    # linear fit on log|dy/dx|-log(x)
    sign = np.sign(np.mean(deriv_data["dydx"]))
    outliers0 = deriv_data[deriv_data["dydx"] == 0.0].index
    r = try_curve_fit(
        lambda x, m, c: m * x + c,
        np.log(deriv_data["x"].drop(outliers0)),
        np.log(abs(deriv_data["dydx"].drop(outliers0))),
        p0=[1.0, 1.0],
        badness_threshold=badness_threshold)
    if r is None:
        return None
    outliers = r["outliers"]
    mc0 = r["p"]
    cov = r["cov"]
    b0 = mc0[0] + 1.0
    a0 = sign * np.exp(mc0[1]) / (mc0[0] + 1.0)
    result["logderiv"] = {
        "loglog_params": mc0,
        "loglog_cov": cov,
        "exponent": b0,
        "coefficient": a0,
        "outliers": list(outliers0) + list(outliers),
    }

    c0, var_c0 = scipy.optimize.curve_fit(
        lambda x, c: a0 * x ** b0 + c,
        data["x"], data["y"], p0=[0.0], maxfev=maxfev)
    result["fixedab"] = {
        "exponent": b0,
        "coefficient": a0,
        "constant": c0[0],
        "constant_err": var_c0[0, 0] ** 0.5,
    }

    try:
        abc, var_abc = scipy.optimize.curve_fit(
            lambda x, a, b, c: a * x ** b + c,
            data["x"], data["y"], p0=[a0, b0, c0], maxfev=maxfev)
    except (RuntimeError, TypeError) as e:
        print(e)
        pass
    else:
        result["full"] = {
            "exponent": abc[1],
            "exponent_err": var_abc[1, 1] ** 0.5,
            "coefficient": abc[0],
            "coefficient_err": var_abc[0, 0] ** 0.5,
            "constant": abc[2],
            "constant_err": var_abc[2, 2] ** 0.5,
            "covariance": var_abc,
        }

    return result

def differentiate(d, xcol, ycol, dydxcol):
    d = d.sort_values(xcol)
    dx = d[xcol].diff()
    dy = d[ycol].diff()
    zero_dx = dx == 0.0
    if np.any(zero_dx):
        if np.all(zero_dx ==
                  (dy.abs() / d[ycol].mean() < DUPLICATE_THRESHOLD)):
            dx = dx[~zero_dx]
            dy = dy[~zero_dx]
        else:
            raise ValueError("data contains irreconcilable duplicates")
    return pd.DataFrame({
        xcol: d[xcol] - dx / 2.0,       # centered x values
        dydxcol: dy / dx,               # dy/dx at centered x values
    }).dropna()

def fit_label(stage, exponent, exponent_err, **kwargs):
    b_err = ""
    if exponent_err is not None:
        b_err = "±{:.2g}".format(exponent_err)
    return "{stage} fit (b={exponent:.3g}{b_err})".format(**locals())

def plot_fits(plot,
              data,
              get_fit_range,
              badness_threshold,
              x_col,
              x_label,
              y_col,
              y_label,
              absdydx_label,
              title_cols,
              get_title,
              get_fn,
              color_cols,
              get_color,
              get_color_label,
              get_dmc,
              dmc_label,
              dmc_yerr_col=None):
    # continuous x range for plotting continuous functions
    x_c = np.linspace(data[x_col].min() - 1, data[x_col].max() + 1, 250)
    [(title_key, gg)] = utils.groupby(data, title_cols)
    fig, ax = plt.subplots(2)
    fig.set_size_inches(8, 10) # otherwise the text will get obscured
    y_range = np.array([np.nan, np.nan])
    fit_results = {}
    for color_key, g in utils.groupby(gg, color_cols):
        logging.info(f"plot_fits: method: {color_key['method']}")
        color = get_color(color_key["method"])
        label = get_color_label(color_key["method"])
        g = g.sort_values([x_col])
        d = g.rename(columns={x_col: "x", y_col: "y"})
        deriv_d = differentiate(d, "x", "y", "dydx")
        x = deriv_d["x"]
        dydx = deriv_d["dydx"]
        ax[0].plot(x, abs(dydx), "x", label=label, color=color)
        ax[1].plot(d["x"], d["y"], "x", label=label, color=color)
        utils.update_range(y_range, d["y"])

        fit_range = get_fit_range(color_key["method"])
        fit_range = (max(fit_range[0], d["x"].min()),
                     min(fit_range[1], d["x"].max()))
        if fit_range[1] < fit_range[0]:
            continue
        fit_range = fit_range + np.array([-0.2, 0.2]) # to make it look nicer
        ax[0].axvspan(fit_range[0], fit_range[1], alpha=0.05, color=color)
        ax[1].axvspan(fit_range[0], fit_range[1], alpha=0.05, color=color)

        d_subset = d[d["x"].between(*fit_range)]
        deriv_subset = deriv_d[deriv_d["x"].between(*fit_range)]
        if len(deriv_subset) < 2:
            continue

        fit = do_fit(d_subset, deriv_subset,
                     badness_threshold=badness_threshold)
        if fit is None:
            continue
        fit_result = {
            "num_points": len(d["x"])
        }
        fit_result.update(fit)
        fit_results[color_key["method"]] = fit_result

        outliers = deriv_subset.loc[fit["logderiv"]["outliers"]]
        ax[0].plot(outliers["x"], abs(outliers["dydx"]), "o",
                   markerfacecolor="none",
                   label="", color="red")
        for stage, result in fit.items():
            if stage == "fixedab":
                continue # fixedab yields the same plot here as logderiv
            a = result["coefficient"]
            b = result["exponent"]
            b_err = result.get("exponent_err", None)
            dydx_c = a * b * x_c ** (b - 1.0)
            ax[0].plot(
                x_c, abs(dydx_c), linestyle=STAGE_TO_LINESTYLE[stage],
                label=label + " " + fit_label(stage, b, b_err),
                color=color)

        for stage, result in fit.items():
            if "constant" not in result:
                continue
            a = result["coefficient"]
            b = result["exponent"]
            c = result["constant"]
            b_err = result.get("exponent_err", None)
            y_c = a * x_c ** b + c
            ax[1].plot(
                x_c, y_c, linestyle=STAGE_TO_LINESTYLE[stage],
                label=label + " " + fit_label(stage, b, b_err),
                color=color)
            if b < 0:
                ax[1].axhline(c, linestyle=":", color=color)
            else:
                logging.warn(f"plot_fits: {stage}.b >= 0: no asymptotic result")
            utils.update_range(y_range, c)

    g = get_dmc(**title_key)
    if len(g):
        y = g[y_col].iloc[0]
        if dmc_yerr_col is not None:
            y_err = g[dmc_yerr_col].iloc[0]
            ax[1].axhspan(y - y_err, y + y_err, alpha=0.4,
                        color="black", label=dmc_label)
            utils.update_range(y_range, [y - y_err, y + y_err])
        # add an extra line to make sure it's visible
        ax[1].axhline(y, alpha=0.4, color="black")
        utils.update_range(y_range, [y])

    ax[0].set_xlabel(x_label)
    ax[0].set_ylabel(absdydx_label)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_title(get_title(**title_key))
    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width * 0.6, box.height])
    ax[0].legend(bbox_to_anchor=(1, 1.0))

    ax[1].legend()
    ax[1].set_xlabel(x_label)
    ax[1].set_ylabel(y_label)
    ax[1].set_ylim(*utils.expand_range(y_range, 0.05))
    box = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0, box.width * 0.6, box.height])
    ax[1].legend(bbox_to_anchor=(1, 1.0))
    ax[1].get_xaxis().set_major_locator(
        matplotlib.ticker.MaxNLocator(integer=True))

    if plot:
        fn = get_fn(**title_key)
        settings_fn = os.path.join("plot_settings", fn + ".json")
        settings = utils.load_json(settings_fn) or {"ax1": {}, "ax2": {}}
        fit_results_fn = os.path.join("fit_results", fn + ".json")
        def save_settings():
            utils.save_json(settings_fn, settings)
        utils.sync_axes_lims(ax[0], settings["ax1"], save_settings)
        utils.sync_axes_lims(ax[1], settings["ax2"], save_settings)
        utils.savefig(fig, fn)
    return fit_results

def plot(label, freq, num_filled, fit_start, fit_stop=np.inf,
         fit_ranges={}, interaction="normal", methods=None,
         badness_threshold=LOGDERIV_BADNESS_THRESHOLD, plot=True):
    '''label: ground, add, or rm.

    fit_ranges: {method: (start, stop)}
    Used to adjust the fit_range of a specific method (overrides fit_start and
    fit_stop).'''

    d_dmc = utils.load_all_dmc()
    d = utils.load_all()
    d = utils.filter_preferred_ml(d)
    d = d[~d["method"].isin(["imsrg[f]+eom[n]"])]

    # filters
    d = d[d["interaction"] == interaction]
    d = d[d["num_filled"] == num_filled]
    d = d[d["freq"] == freq]
    d = d[d["label"] == label]
    if methods is not None:
        d = d[d["method"].isin(methods)]

    if label == "ground":
        e_sym = "E/N"
        e_text = "energy per particle"
        d["energy"] /= d["num_particles"]
        d_dmc["energy"] /= d_dmc["num_particles"]
        d_dmc["energy_err"] /= d_dmc["num_particles"]
    else:
        e_sym = "ε"
        e_text = "energy"

    default_fit_range = (fit_start, fit_stop)
    fit_results = plot_fits(
        plot=plot,
        data=d,
        get_fit_range=lambda method: fit_ranges.get(method, default_fit_range),
        badness_threshold=badness_threshold,
        x_col="num_shells",
        x_label="K (number of shells)",
        y_col="energy",
        y_label=f"{e_sym} ({e_text})",
        absdydx_label=f"|∂{e_sym}/∂K|",
        title_cols=["label", "num_filled", "num_particles", "freq"],
        get_title="{label}, N={num_particles}, ω={freq}".format,
        get_fn="fit-{label}-{num_filled}-{freq}".format,
        color_cols=["method"],
        get_color=utils.METHOD_COLOR.get,
        get_color_label=utils.METHOD_LABEL.get,
        get_dmc=utils.filter_eq(d_dmc, ["label", "num_filled", "freq"],
                                check_unused_kwargs=False),
        dmc_label="DMC",
        dmc_yerr_col="energy_err",
    )
    results = []
    for method, fit_result in sorted(fit_results.items()):
        results.append({
            "interaction": interaction,
            "label": label,
            "freq": freq,
            "num_filled": num_filled,
            "method": method,
            "fit_result": fit_result,
        })
    return results

def main():
    r = []

    r.extend(plot("ground", num_filled=1, freq=1.0, fit_start=6))
    r.extend(plot("ground", num_filled=1, freq=0.5, fit_start=6))
    r.extend(plot("ground", num_filled=1, freq=0.28, fit_start=6))
    r.extend(plot("ground", num_filled=1, freq=0.1, fit_start=8))
    r.extend(plot("ground", num_filled=2, freq=1.0, fit_start=6))
    r.extend(plot("ground", num_filled=2, freq=0.5, fit_start=6))
    r.extend(plot("ground", num_filled=2, freq=0.28, fit_start=6))
    r.extend(plot("ground", num_filled=2, freq=0.1, fit_start=10))
    r.extend(plot("ground", num_filled=3, freq=1.0, fit_start=9))
    r.extend(plot("ground", num_filled=3, freq=0.5, fit_start=9))
    r.extend(plot("ground", num_filled=3, freq=0.28, fit_start=11))
    r.extend(plot("ground", num_filled=3, freq=0.1, fit_start=11))
    r.extend(plot("ground", num_filled=4, freq=1.0, fit_start=10))
    r.extend(plot("ground", num_filled=4, freq=0.5, fit_start=12,
                  fit_ranges={"hf": [10, 99999]}))
    r.extend(plot("ground", num_filled=4, freq=0.28, fit_start=14))
    r.extend(plot("ground", num_filled=4, freq=0.1, fit_start=11))
    r.extend(plot("ground", num_filled=5, freq=1.0, fit_start=13))
    r.extend(plot("ground", num_filled=5, freq=0.5, fit_start=14))
    r.extend(plot("ground", num_filled=5, freq=0.28, fit_start=17))
    r.extend(plot("ground", num_filled=5, freq=0.1, fit_start=11))
    r.extend(plot("ground", num_filled=6, freq=1.0, fit_start=16))
    r.extend(plot("ground", num_filled=6, freq=0.28, fit_start=15))
    r.extend(plot("ground", num_filled=6, freq=0.1, fit_start=15,
                  badness_threshold=0.0))

    r.extend(plot("add", num_filled=1, freq=1.0, fit_start=7))
    r.extend(plot("add", num_filled=1, freq=0.5, fit_start=7))
    r.extend(plot("add", num_filled=1, freq=0.28, fit_start=7))
    r.extend(plot("add", num_filled=2, freq=1.0, fit_start=9))
    r.extend(plot("add", num_filled=2, freq=0.5, fit_start=9))
    r.extend(plot("add", num_filled=2, freq=0.28, fit_start=10))
    r.extend(plot("add", num_filled=2, freq=0.1, fit_start=10))
    r.extend(plot("add", num_filled=3, freq=1.0, fit_start=13))
    r.extend(plot("add", num_filled=3, freq=0.5, fit_start=10))
    r.extend(plot("add", num_filled=3, freq=0.28, fit_start=13))
    r.extend(plot("add", num_filled=3, freq=0.1, fit_start=10))
    r.extend(plot("add", num_filled=4, freq=1.0, fit_start=13))
    r.extend(plot("add", num_filled=4, freq=0.28, fit_start=13))
    r.extend(plot("add", num_filled=4, freq=0.1, fit_start=11))
    r.extend(plot("add", num_filled=5, freq=1.0, fit_start=10))
    r.extend(plot("add", num_filled=5, freq=0.5, fit_start=10))
    r.extend(plot("add", num_filled=5, freq=0.28, fit_start=10))
    r.extend(plot("add", num_filled=5, freq=0.1, fit_start=10))

    r.extend(plot("rm", num_filled=1, freq=1.0, fit_start=7))
    r.extend(plot("rm", num_filled=1, freq=0.5, fit_start=7))
    r.extend(plot("rm", num_filled=1, freq=0.28, fit_start=7,
                  fit_ranges={"hf+qdpt3": [12, 99999]}))
    r.extend(plot("rm", num_filled=2, freq=1.0, fit_start=10))
    r.extend(plot("rm", num_filled=2, freq=0.5, fit_start=7))
    r.extend(plot("rm", num_filled=2, freq=0.28, fit_start=10))
    r.extend(plot("rm", num_filled=2, freq=0.1, fit_start=9))
    r.extend(plot("rm", num_filled=3, freq=1.0, fit_start=11))
    r.extend(plot("rm", num_filled=3, freq=0.5, fit_start=9))
    r.extend(plot("rm", num_filled=3, freq=0.28, fit_start=11))
    r.extend(plot("rm", num_filled=3, freq=0.1, fit_start=9))
    r.extend(plot("rm", num_filled=4, freq=1.0, fit_start=13))
    r.extend(plot("rm", num_filled=4, freq=0.5, fit_start=10))
    r.extend(plot("rm", num_filled=4, freq=0.28, fit_start=10))
    r.extend(plot("rm", num_filled=4, freq=0.1, fit_start=10))
    r.extend(plot("rm", num_filled=5, freq=1.0, fit_start=13))
    r.extend(plot("rm", num_filled=5, freq=0.5, fit_start=10))
    r.extend(plot("rm", num_filled=5, freq=0.28, fit_start=10))
    r.extend(plot("rm", num_filled=5, freq=0.1, fit_start=11))

    utils.save_json("fit_results.json", r)

if __name__ == "__main__":
    warnings.simplefilter("ignore", scipy.optimize.OptimizeWarning)
    utils.plot_main(__file__, plot, main)
