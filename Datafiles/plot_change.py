#!/usr/bin/env python3
# Generates ../FigureFiles/fig-change-*.svg from the data files.
import argparse, json, os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize
import utils

FIT_TYPES = ["loglog", "semilog"]

def plot(interactive, fit_type):
    if interactive:
        plt.ion()

    suffix = ""
    fit_points = 4

    if fit_type not in FIT_TYPES:
        raise ValueError("fit_type must be one of: {}".format(FIT_TYPES))

    d = pd.concat(utils.get_ar_energies_for_v(suffix), ignore_index=True)

    d = d[d["num_filled"].isin([2])]
    d = d[d["freq"].isin([0.28, 1.0])]
    d = d[((d["num_shells"] >= 4) &
           (d["num_shells"] <= 15))]
    d = d[d["method"].isin(["qdpt", "eom", "eom_quads", "cc"])]

    # use the ml closest to shell closure
    d = d[d["ml"] == (d["num_filled"] + (d["label"] != "add")) % 2]

    d["num_particles"] = d["num_filled"] * (d["num_filled"] + 1)

    d["num_orbitals"] = d["num_shells"] * (d["num_shells"] + 1)

    d["title"] = d.apply(lambda r: "{label}, n={num_particles}, ω={freq}"
                         .format(**r), axis=1)
    d["filename"] = d.apply(lambda r: os.path.join(
        os.path.pardir,
        "FigureFiles",
        "fig-change-{num_particles}-{freq}-{label}-{ml}.svg"
        .format(**r),
    ), axis=1)

    d["x"] = d["num_orbitals"]
    d["y"] = d["energy"]

    # calculating finite differences assumes a certain ordering
    d = d.sort_values(["x"])
    gs = d.groupby(["title", "method"])
    dy = pd.concat(g["y"].diff() for _, g in gs)
    dx = pd.concat(g["x"].diff() for _, g in gs)

    # replace y with its approximate derivative and shift x appropriately
    d["y"] = abs(dy / dx)
    d["x"] = d["x"] - dx / 2.0
    d = d.dropna()

    def fit(x, y, **params):
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
        badness = (np.sum((f(x_in, *p) - y_in) ** 2) / (len(x_in) - len(p)))
        r = {
            "fit_type": fit_type,
            "badness": float(badness),
            "params": p.tolist(),
            "cov_params": var_p.tolist(),
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
        sys.stdout.write("{}\n\n".format(
            json.dumps(r, **utils.JSON_PRETTY)))
        sys.stdout.flush()
        return x, f(x, *p)

    color_map = {
        "qdpt": "#39b237",
        "eom": "#841a0b",
        "eom_quads": "#a825bc",
        "cc": "#1351c4",
    }

    for title, fig_data in d.groupby("title"):
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 10)

        for method, g in fig_data.groupby("method"):
            ax.plot(g["x"], g["y"], "x",
                    label=method,
                    color=color_map[method])
            fit_x, fit_y = fit(x=g["x"], y=g["y"],
                               title=title, method=method)
            ax.plot(fit_x, fit_y, "-",
                    label=(method + " fit"),
                    linewidth=1.5,
                    color=color_map[method])
        ax.legend()
        ax.set_xlabel("N (number of orbitals)")
        ax.set_xscale("log" if fit_type == "loglog" else "linear")
        ax.set_ylabel("(dε/dN)/ε")
        ax.set_yscale("log")
        ax.set_title(title)

        if not interactive:
            # FIXME: not elegant
            fn, = fig_data["filename"].unique()
            fig.savefig(fn)
            plt.close(fig)
            sys.stderr.write("// Figure saved to: {}\n\n".format(fn))
            sys.stderr.flush()

    if interactive:
        plt.show(block=True)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--interactive", action="store_true")
    kwargs = vars(p.parse_args())
    utils.init(__file__)
    for fit_type in FIT_TYPES:
        plot(fit_type=fit_type, **kwargs)

main()
