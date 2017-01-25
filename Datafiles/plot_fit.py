#!/usr/bin/env python3
import argparse, itertools, json, os, scipy, sys, traceback, warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utils

DUPLICATE_THRESHOLD = 1e-6

STAGE_TO_LINESTYLE = {
    "logderiv": ":",
    "fixedab": "--",
    "full": "-",
}

def default_fit_condition(cov, **kwargs):
    return np.isfinite(cov).all()

def default_fit_compare(cov1, cov2, **kwargs):
    # dimensionally this is valid and has the advantage of not being dependent
    # on the parameter values, but I'm not sure if there is a rigorous
    # justification to its sensibleness
    return -np.sign(np.linalg.norm(cov1 / cov2) - cov1.shape[0] ** 2)

class FittingError(Exception):
    pass

def try_curve_fit(f, x, y, p0,
                  max_num_outliers=2,
                  fit_compare=default_fit_compare):
    '''Like scipy.optimize.curve_fit, but tries to drop some of the points to
    find the most satisfactory fit.  The quality of the fit is compared using
    'fit_compare'.'''
    max_num_outliers = min(max_num_outliers, len(x) - len(p0))
    best_outliers = None
    for num_outliers in range(max_num_outliers + 1):
        for outliers in itertools.combinations(list(x.index), num_outliers):
            outliers = list(outliers)
            new_x = x.drop(outliers)
            new_y = y.drop(outliers)
            p, cov = scipy.optimize.curve_fit(f, new_x, new_y, p0=p0)
            if (best_outliers is None or
                    fit_compare(p1=best_p, cov1=best_cov,
                                outliers1=best_outliers,
                                p2=p, cov2=cov, outliers2=outliers) < 0):
                best_outliers = outliers
                best_p = p
                best_cov = cov
    if best_outliers is not None:
        return best_outliers, best_p, best_cov
    raise FittingError("unable to get a satisfying fit")

def do_fit(data, deriv_data):
    '''Perform a fit using a power law model:

        y = a * x ** b + c

    data: DataFrameish{"x", "y"}
    deriv_data: DataFrameish{"x", "dydx"}
    '''

    result = {}

    # linear fit on log|dy/dx|-log(x)
    sign = np.sign(np.mean(deriv_data["dydx"]))
    outliers0 = deriv_data[deriv_data["dydx"] == 0.0].index
    outliers, mc0, cov = try_curve_fit(
        lambda x, m, c: m * x + c,
        np.log(deriv_data["x"].drop(outliers0)),
        np.log(abs(deriv_data["dydx"].drop(outliers0))),
        p0=[1.0, 1.0])
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
        data["x"], data["y"], p0=[0.0])
    result["fixedab"] = {
        "exponent": b0,
        "coefficient": a0,
        "constant": c0[0],
        "constant_err": var_c0[0] ** 0.5,
    }

    try:
        abc, var_abc = scipy.optimize.curve_fit(
            lambda x, a, b, c: a * x ** b + c,
            data["x"], data["y"], p0=[a0, b0, c0])
    except RuntimeError:
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

def plot_fits(data,
              fit_range,
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
    if fit_range[1] == np.inf:
        fit_range[1] = np.max(data[x_col])
    fit_range = fit_range + np.array([-0.2, 0.2]) # make it look less ambiguous
    # continuous x range for plotting continuous functions
    x_c = np.linspace(data[x_col].min() - 1, data[x_col].max() + 1, 250)
    [(title_key, gg)] = utils.groupby(data, title_cols)
    fig, ax = plt.subplots(2)
    fig.set_size_inches(8, 10) # otherwise the text will get obscured
    y_range = np.array([np.nan, np.nan])
    fit_results = {}
    for color_key, g in utils.groupby(gg, color_cols):
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

        d_subset = d[d["x"].between(*fit_range)]
        deriv_subset = deriv_d[deriv_d["x"].between(*fit_range)]
        if len(deriv_subset) < 2:
            continue
        try:
            fit = do_fit(d_subset, deriv_subset)
        except FittingError as e:
            sys.stderr.write(f"could not fit {color_key['method']}: {e}\n")
            sys.stderr.flush()
            continue

        fit_results[color_key["method"]] = fit

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
            ax[1].axhline(c, linestyle=":", color=color)
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

    ax[0].axvspan(fit_range[0], fit_range[1], alpha=0.15, color="#d6a528")
    ax[0].set_xlabel(x_label)
    ax[0].set_ylabel(absdydx_label)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_title(get_title(**title_key))
    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width * 0.6, box.height])
    ax[0].legend(bbox_to_anchor=(1, 1.0))

    ax[1].axvspan(fit_range[0], fit_range[1], alpha=0.15, color="#d6a528")
    ax[1].legend()
    ax[1].set_xlabel(x_label)
    ax[1].set_ylabel(y_label)
    ax[1].set_ylim(*utils.expand_range(y_range, 0.05))
    box = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0, box.width * 0.6, box.height])
    ax[1].legend(bbox_to_anchor=(1, 1.0))

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
         interaction="normal", hartree_fock=False):
    d_dmc = utils.load_all_dmc()
    d = utils.load_all()
    d = utils.filter_preferred_ml(d)
    d = d[d["method"] != "imsrg[f]+eom[n]"]

    # filters
    d = d[d["interaction"] == interaction]
    d = d[d["num_filled"] == num_filled]
    d = d[d["freq"] == freq]
    d = d[d["label"] == label]
    if hartree_fock:
        d = d[d["method"] != "hf"]

    if label == "ground":
        e_sym = "E/N"
        e_text = "energy per particle"
        d["energy"] /= d["num_particles"]
        d_dmc["energy"] /= d_dmc["num_particles"]
        d_dmc["energy_err"] /= d_dmc["num_particles"]
    else:
        e_sym = "ε"
        e_text = "energy"
    fit_results = plot_fits(
        data=d,
        fit_range=[fit_start, fit_stop],
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
    r.extend(plot("ground", num_filled=1, freq=0.28, fit_start=6))
    r.extend(plot("ground", num_filled=2, freq=1.0, fit_start=6))
    r.extend(plot("ground", num_filled=2, freq=0.28, fit_start=6))
    r.extend(plot("ground", num_filled=3, freq=1.0, fit_start=9))
    r.extend(plot("ground", num_filled=3, freq=0.28, fit_start=9))
    r.extend(plot("ground", num_filled=4, freq=1.0, fit_start=10))
    r.extend(plot("ground", num_filled=4, freq=0.28, fit_start=14))
    r.extend(plot("ground", num_filled=5, freq=1.0, fit_start=15))
    r.extend(plot("ground", num_filled=5, freq=0.28, fit_start=17))
    r.extend(plot("ground", num_filled=6, freq=1.0, fit_start=16))

    r.extend(plot("add", num_filled=1, freq=1.0, fit_start=7))
    r.extend(plot("add", num_filled=1, freq=0.28, fit_start=7))
    r.extend(plot("add", num_filled=2, freq=1.0, fit_start=10))
    r.extend(plot("add", num_filled=2, freq=0.28, fit_start=10))
    r.extend(plot("add", num_filled=3, freq=1.0, fit_start=13))
    r.extend(plot("add", num_filled=3, freq=0.28, fit_start=13))
    r.extend(plot("add", num_filled=4, freq=1.0, fit_start=13))

    r.extend(plot("rm", num_filled=1, freq=1.0, fit_start=7))
    r.extend(plot("rm", num_filled=1, freq=0.28, fit_start=7))
    r.extend(plot("rm", num_filled=2, freq=1.0, fit_start=10))
    r.extend(plot("rm", num_filled=2, freq=0.28, fit_start=10))
    r.extend(plot("rm", num_filled=3, freq=1.0, fit_start=13))
    r.extend(plot("rm", num_filled=3, freq=0.28, fit_start=13))
    r.extend(plot("rm", num_filled=4, freq=1.0, fit_start=13))

    utils.save_json("fit_results.json", r)

warnings.simplefilter("ignore", scipy.optimize.OptimizeWarning)
utils.plot_main(__file__, plot, main)
