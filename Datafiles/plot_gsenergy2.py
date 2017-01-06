#!/usr/bin/env python3
import argparse, json, os, scipy, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utils

def do_fit(data, deriv_data):
    '''Perform a fit using a power law model

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
    result["logderiv"] = r["extra0"] # linear fit on log|dy/dx|-log(x)
    try:
        result["deriv"] = r["extra"]     # power fit on dy/dx-x
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

def plot(interactive):

    d = pd.read_csv("imsrg-qdpt/dat_gsenergy.txt",
                    header=0, index_col=False,
                    delim_whitespace=True)

    d_dmc = pd.read_csv("dat-qdsfe.jEYRh-4ptC1Dr5nlBcz0tg/dat-gs-dmc.txt",
                        header=0, index_col=False,
                        delim_whitespace=True, comment="#")

    d = d[d["num_filled"] == 6]
    d = d[d["freq"] == 1.0]
    d = d[d["method"] != "hf"]
    fit_range = [16, 20]

    # TODO: make use of this
    preferred_fit_ranges = [
        ["num_filled", "freq", "fit_range_start"],
        [1, 1.0, 6],
        [1, 0.28, 6],
        [2, 1.0, 6],
        [2, 0.28, 6],
        [3, 1.0, 9],
        [3, 0.28, 9],
        [4, 1.0, 10],
        # 4, 0.28, ??? # TODO: need more data
        [5, 1.0, 15],
        [5, 0.28, 17],
        [6, 1.0, 16],
    ]
    # TODO: seems like IMSRG is better than MP2 for these cases
    #       but for larger num_filled >= 6 cases it seems worse?
    #       need more data!

    linestyle = {
        "logderiv": ":",
        "deriv": "-.",
        "fixedab": "--",
        "full": "-",
    }

    fit_range = fit_range + np.array([-0.2, 0.2]) # make it look less ambiguous
    for [num_filled, freq], gg in d.groupby(["num_filled", "freq"]):
        num_particles = num_filled * (num_filled + 1)
        fig1, ax1 = plt.subplots()
        fig1.set_size_inches(8, 10)
        fig2, ax2 = plt.subplots()
        fig2.set_size_inches(8, 10)
        y_range = np.array([np.nan, np.nan])

        for method, g in gg.groupby(["method"]):
            color = utils.GS_METHOD_COLOR[method]

            g = g.sort_values(["num_shells"])
            data = g.rename(columns={
                "num_shells": "x",
                "energy": "y",
            })
            dx = data["x"].diff()
            dy = data["y"].diff()
            deriv_data = pd.DataFrame({
                "x": data["x"] - dx / 2.0, # centered x values
                "dydx": dy / dx,  # dy/dx at centered x values
            }).dropna()

            fit = do_fit(
                data[data["x"].between(*fit_range)],
                deriv_data[deriv_data["x"].between(*fit_range)])
            sys.stdout.write(utils.json_pretty(fit) + "\n\n")
            sys.stdout.flush()

            def fit_label(method, stage, exponent, exponent_err, **kwargs):
                b_err = ""
                if exponent_err is not None:
                    b_err = "±{exponent_err:.2g}".format(**locals())
                return ("{method} {stage} fit (b={exponent:.3g}{b_err})"
                        .format(**locals()))

            x = deriv_data["x"]
            dydx = deriv_data["dydx"]
            ax1.plot(x, abs(dydx), "o", label=method, color=color)
            for stage, result in fit.items():
                if stage == "fixedab":
                    continue # fixedab yields the same plot here as logderiv
                a = result["coefficient"]
                b = result["exponent"]
                b_err = result.get("exponent_err", None)
                dydx = a * b * x ** (b - 1.0)
                ax1.plot(x, abs(dydx), linestyle=linestyle[stage],
                        label=fit_label(method, stage, b, b_err),
                        color=color)

            x = data["x"]
            y = data["y"]
            ax2.plot(x, y / num_particles, "o", label=method, color=color)
            utils.update_range(y_range, y)
            for stage, result in fit.items():
                if "constant" not in result:
                    continue
                a = result["coefficient"]
                b = result["exponent"]
                c = result["constant"]
                b_err = result.get("exponent_err", None)
                y = a * x ** b + c
                ax2.plot(x, y / num_particles,
                         linestyle=linestyle[stage],
                        label=fit_label(method, stage, b, b_err),
                        color=color)
                ax2.axhline(c / num_particles, linestyle=":", color=color)
                utils.update_range(y_range, c)

        g = d_dmc[(d_dmc["num_filled"] == num_filled) &
                  d_dmc["freq"].between(freq - 1e-8, freq + 1e-8)]
        if len(g):
            y = g["energy"].iloc[0]
            y_err = g["energy_err"].iloc[0]
            ax2.axhspan((y - y_err) / num_particles,
                        (y + y_err) / num_particles,
                        alpha=0.4,
                        color="black",
                        label="DMC")
            # add an extra line to make sure it's visible
            ax2.axhline(y / num_particles, alpha=0.4, color="black")
            utils.update_range(y_range, [y - y_err, y + y_err])

        ax1.axvspan(fit_range[0], fit_range[1], alpha=0.15, color="#d6a528")
        ax1.legend()
        ax1.set_xlabel("K (number of shells)")
        ax1.set_ylabel("|ΔE/ΔK|")
        ax1.set_xscale("log")
        ax1.set_yscale("log")

        ax2.axvspan(fit_range[0], fit_range[1], alpha=0.15, color="#d6a528")
        ax2.legend()
        ax2.set_xlabel("K (number of shells)")
        ax2.set_ylabel("E/N (energy per particle)")

        ax2.set_ylim(*utils.expand_range(y_range / num_particles, 0.05))

        utils.savefig(
            fig1, interactive,
            "gsenergy2-change-{num_filled}-{freq}".format(**locals()))
        utils.savefig(
            fig2, interactive,
            "gsenergy2-{num_filled}-{freq}".format(**locals()))

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--interactive", action="store_true")
    kwargs = vars(p.parse_args())
    with utils.plot(__file__, interactive=kwargs["interactive"]):
        plot(**kwargs)

main()
