#!/usr/bin/env python3
import argparse, json, os, scipy, sys, traceback
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utils

DUPLICATE_THRESHOLD = 1e-6

STAGE_TO_LINESTYLE = {
    "logderiv": ":",
    "deriv": "-.",
    "fixedab": "--",
    "full": "-",
}

def do_fit(data, deriv_data):
    '''Perform a fit using a power law model:

        y = a * x ** b + c

    data: DataFrameish{"x", "y"}
    deriv_data: DataFrameish{"x", "dydx"}
    '''
    result = {}

    r = utils.fit_change(
        fit_type="loglog",
        fit_points=0,
        x=deriv_data["x"],
        y=deriv_data["dydx"],
    )
    result["logderiv"] = r["extra0"]    # linear fit on log|dy/dx|-log(x)
    try:
        result["deriv"] = r["extra"]    # power fit on dy/dx-x
    except KeyError:
        pass

    b0 = result["logderiv"]["exponent"]
    a0 = result["logderiv"]["coefficient"]
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
    x_c = np.linspace(data[x_col].min(), data[x_col].max(), 250) # continuous x
    for title_key, gg in utils.groupby(data, title_cols):
        fig, [ax1, ax2] = plt.subplots(2)
        fig.set_size_inches(8, 10) # otherwise the text will get obscured
        y_range = np.array([np.nan, np.nan])
        for color_key, g in utils.groupby(gg, color_cols):
            color = get_color(**color_key)
            label = get_color_label(**color_key)
            g = g.sort_values([x_col])
            d = g.rename(columns={x_col: "x", y_col: "y"})
            deriv_d = differentiate(d, "x", "y", "dydx")
            x = deriv_d["x"]
            dydx = deriv_d["dydx"]
            ax1.plot(x, abs(dydx), "x", label=label, color=color)

            d_subset = d[d["x"].between(*fit_range)]
            deriv_subset = deriv_d[deriv_d["x"].between(*fit_range)]
            if len(deriv_subset) < 2:
                continue
            fit = do_fit(d_subset, deriv_subset)
            sys.stdout.write(utils.json_pretty(fit) + "\n\n")
            sys.stdout.flush()

            for stage, result in fit.items():
                if stage == "fixedab":
                    continue # fixedab yields the same plot here as logderiv
                a = result["coefficient"]
                b = result["exponent"]
                b_err = result.get("exponent_err", None)
                dydx_c = a * b * x_c ** (b - 1.0)
                ax1.plot(
                    x_c, abs(dydx_c), linestyle=STAGE_TO_LINESTYLE[stage],
                    label=label + " " + fit_label(stage, b, b_err),
                    color=color)

            x = d["x"]
            y = d["y"]
            ax2.plot(x, y, "x", label=label, color=color)
            utils.update_range(y_range, y)
            for stage, result in fit.items():
                if "constant" not in result:
                    continue
                a = result["coefficient"]
                b = result["exponent"]
                c = result["constant"]
                b_err = result.get("exponent_err", None)
                y_c = a * x_c ** b + c
                ax2.plot(
                    x_c, y_c, linestyle=STAGE_TO_LINESTYLE[stage],
                    label=label + " " + fit_label(stage, b, b_err),
                    color=color)
                ax2.axhline(c, linestyle=":", color=color)
                utils.update_range(y_range, c)

        g = get_dmc(**title_key)
        if len(g):
            y = g[y_col].iloc[0]
            if dmc_yerr_col is not None:
                y_err = g[dmc_yerr_col].iloc[0]
                ax2.axhspan(y - y_err, y + y_err, alpha=0.4,
                            color="black", label=dmc_label)
                utils.update_range(y_range, [y - y_err, y + y_err])
            # add an extra line to make sure it's visible
            ax2.axhline(y, alpha=0.4, color="black")
            utils.update_range(y_range, [y])

        ax1.axvspan(fit_range[0], fit_range[1], alpha=0.15, color="#d6a528")
        ax1.set_xlabel(x_label)
        ax1.set_ylabel(absdydx_label)
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_title(get_title(**title_key))
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.6, box.height])
        ax1.legend(bbox_to_anchor=(1, 1.0))

        ax2.axvspan(fit_range[0], fit_range[1], alpha=0.15, color="#d6a528")
        ax2.legend()
        ax2.set_xlabel(x_label)
        ax2.set_ylabel(y_label)
        ax2.set_ylim(*utils.expand_range(y_range, 0.05))
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.6, box.height])
        ax2.legend(bbox_to_anchor=(1, 1.0))

        fn = get_fn(**title_key)
        settings_fn = fn + ".settings"
        settings = utils.load_json(settings_fn) or {"ax1": {}, "ax2": {}}
        def save_settings():
            utils.save_json(settings_fn, settings)
        utils.sync_axes_lims(ax1, settings["ax1"], save_settings)
        utils.sync_axes_lims(ax2, settings["ax2"], save_settings)
        utils.savefig(fig, fn)
        break

def plot_addrm(label, num_filled, freq, fit_start, fit_stop, **kwargs):
    d_dmc = pd.read_csv("dat-qdsfe.jEYRh-4ptC1Dr5nlBcz0tg/dat-addrm-dmc.txt",
                        header=0, index_col=False,
                        delim_whitespace=True, comment="#")
    d = pd.concat(utils.get_ar_energies_for_v(v=""), ignore_index=True)
    # d = d[d["method"] != "eom_f"]
    # use the ml closest to shell closure
    d = d[d["ml"] == (d["num_filled"] + (d["label"] != "add")) % 2]
    d["num_particles"] = d["num_filled"] * (d["num_filled"] + 1)

    # filters
    d = d[d["num_filled"] == num_filled]
    d = d[d["freq"] == freq]
    d = d[d["label"] == label]

    plot_fits(
        data=d,
        fit_range=[fit_start, fit_stop],
        x_col="num_shells",
        x_label="K (number of shells)",
        y_col="energy",
        y_label="E (energy)",
        absdydx_label="|∂E/∂K|",
        title_cols=["label", "num_filled", "num_particles", "freq"],
        get_title="{label}, n={num_particles}, ω={freq}".format,
        get_fn="fit-{label}-{num_filled}-{freq}".format,
        color_cols=["method"],
        get_color=lambda method: utils.METHOD_COLOR[method],
        get_color_label=lambda method: method,
        get_dmc=lambda **k: d_dmc[
            (d_dmc["is_hole"] == (k["label"] == 0)) &
            (d_dmc["num_filled"] == k["num_filled"]) &
            (d_dmc["freq"] == k["freq"])],
        dmc_label="DMC",
    )

def plot_ground(num_filled, freq, fit_start, fit_stop, **kwargs):
    d_dmc = pd.read_csv("dat-qdsfe.jEYRh-4ptC1Dr5nlBcz0tg/dat-gs-dmc.txt",
                        header=0, index_col=False,
                        delim_whitespace=True, comment="#")
    d = pd.read_csv("imsrg-qdpt/dat_gsenergy.txt",
                    header=0, index_col=False,
                    delim_whitespace=True)
    for data in [d, d_dmc]:
        data["num_particles"] = data["num_filled"] * (data["num_filled"] + 1)
        data["energy_per_particle"] = data["energy"] / data["num_particles"]
    d_dmc["energy_per_particle_err"] = (d_dmc["energy_err"] /
                                        d_dmc["num_particles"])

    # filters
    d = d[d["num_filled"] == num_filled]
    d = d[d["freq"] == freq]
    d = d[d["method"] != "hf"]

    plot_fits(
        data=d,
        fit_range=[fit_start, fit_stop],
        x_col="num_shells",
        x_label="K (number of shells)",
        y_col="energy_per_particle",
        y_label="E/n (energy per particle)",
        absdydx_label="|∂(E/n)/∂K|",
        title_cols=["num_filled", "num_particles", "freq"],
        get_title="ground, n={num_particles}, ω={freq}".format,
        get_fn="fit-ground-{num_filled}-{freq}".format,
        color_cols=["method"],
        get_color=lambda method: utils.GS_METHOD_COLOR[method],
        get_color_label=lambda method: method,
        get_dmc=lambda **k: d_dmc[
            (d_dmc["num_filled"] == k["num_filled"]) &
            (d_dmc["freq"] == k["freq"])],
        dmc_label="DMC/exact",
        dmc_yerr_col="energy_per_particle_err",
    )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--interactive", action="store_true")
    p.add_argument("-fq", "--freq", type=float, required=True)
    p.add_argument("-kf", "--num-filled", type=int, required=True)
    p.add_argument("-k1", "--fit-start", type=float, required=True)
    p.add_argument("-k2", "--fit-stop", type=float, default=np.inf)
    p.add_argument("label", metavar="type", help="ground, add, or rm")
    kwargs = vars(p.parse_args())
    plt.rcParams["interactive"] = kwargs["interactive"]
    with utils.plot(__file__):
        if kwargs["label"] == "ground":
            plot_ground(**kwargs)
        elif kwargs["label"] in ["add", "rm"]:
            plot_addrm(**kwargs)
        else:
            raise ValueError("<type> argument is invalid")

main()
