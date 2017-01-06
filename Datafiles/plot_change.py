#!/usr/bin/env python3
# Generates ../FigureFiles/fig-change-*.svg from the data files.
import argparse, json, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import utils

FIT_TYPES = ["loglog", "semilog"]

def plot(interactive, fit_type):
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

    d["title"] = d.apply(lambda r: "{label}, n={num_particles}, ω={freq}"
                         .format(**r), axis=1)
    d["filename"] = d.apply(lambda r: os.path.join(
        os.path.pardir,
        "FigureFiles",
        "fig-change-{num_particles}-{freq}-{label}-{ml}-{fit_type}.svg"
        .format(fit_type=fit_type, **r),
    ), axis=1)

    d["x"] = d["num_shells"]
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

    def fit(**kwargs):
        r = utils.fit_change(fit_type=fit_type,
                             fit_points=fit_points,
                             **kwargs)
        sys.stdout.write("{}\n\n".format(
            json.dumps(utils.sanitize_json(r["extra"]), **utils.JSON_PRETTY)))
        sys.stdout.flush()
        return r["x"], r["y"]

    for title, fig_data in d.groupby("title"):
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 10)

        for method, g in fig_data.groupby("method"):
            ax.plot(g["x"], g["y"], "x",
                    label=method,
                    color=utils.METHOD_COLOR[method])
            fit_x, fit_y = fit(x=g["x"], y=g["y"],
                               title=title, method=method)
            ax.plot(fit_x, fit_y, "-",
                    label=(method + " fit"),
                    linewidth=1.5,
                    color=utils.METHOD_COLOR[method])
        ax.legend()
        ax.set_xlabel("K (number of shells)")
        ax.set_xscale("log" if fit_type == "loglog" else "linear")
        ax.set_ylabel("Δε/ΔK")
        ax.set_yscale("log")
        ax.set_title(title)

        if not interactive:
            # FIXME: not elegant
            fn, = fig_data["filename"].unique()
            fig.savefig(fn)
            plt.close(fig)
            sys.stderr.write("// Figure saved to: {}\n\n".format(fn))
            sys.stderr.flush()

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--interactive", action="store_true")
    kwargs = vars(p.parse_args())
    with utils.plot(__file__, interactive=kwargs["interactive"]):
        for fit_type in FIT_TYPES:
            plot(fit_type=fit_type, **kwargs)

main()
