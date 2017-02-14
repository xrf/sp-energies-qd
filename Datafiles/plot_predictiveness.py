#!/usr/bin/env python3
import collections, functools, itertools, logging, multiprocessing, random
import matplotlib.cm
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize
import scipy.special
import plot_fit, utils
import plot_fit

CMAP_VIRIDIS = matplotlib.cm.get_cmap("viridis")

CMAP_FOLDED_VIRIDIS = utils.colormap(
    "folded_viridis",
    np.concatenate([CMAP_VIRIDIS(np.linspace(0.0, 1.0, 100)),
                    CMAP_VIRIDIS(np.linspace(1.0, 0.0, 100))]))

CMAP_FOLDED_PINK = utils.colormap(
    "folded_pink",
    ["#ffffff", "#ec7ca3", "#ffffff"])

# scaling transformations
TRANSFORM_LOG_ABS = lambda x: np.log(abs(x)), lambda x: np.exp(x)
TRANSFORM_ID = lambda x: x, lambda x: x

def gather_fit_data_inner(group, fit_count, maxfev):
    badness_threshold = plot_fit.LOGDERIV_BADNESS_THRESHOLD
    (label, interaction, num_filled, freq, method), gg = group
    logging.info(f"{(label, interaction, num_filled, freq, method)}")
    gg = gg.rename(columns={"num_shells": "x", "energy": "y"})
    gg = gg.sort_values(["x"])
    deriv_gg = plot_fit.differentiate(gg, "x", "y", "dydx")
    subresults = []
    for fit_start in range(4, 20):
        fit_stop = fit_start + fit_count - 1
        g = gg[gg["x"].between(fit_start, fit_stop)]
        if len(g) != fit_count:
            continue
        deriv_g = deriv_gg[deriv_gg["x"].between(fit_start, fit_stop)]
        fit = plot_fit.do_fit(g, deriv_g, badness_threshold=badness_threshold,
                              maxfev=maxfev)
        if fit is None:
            continue
        main_fit = fit.get("full", fit["fixedab"])
        exponent = main_fit["exponent"]
        if exponent > 0:
            continue
        constant = main_fit["constant"]
        constant_err = main_fit["constant_err"]
        subresults.append({
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
        })
    results = []
    for r in subresults[:-1]:
        best_constant = subresults[-1]["constant"]
        best_constant_err = subresults[-1]["constant_err"]
        constant = r["constant"]
        constant_err = r["constant_err"]
        dist = r["dist"]
        results.append({
            "label": label,
            "interaction": interaction,
            "num_filled": num_filled,
            "freq": freq,
            "method": method,
            "fit_method": r["fit_method"],
            "fit_stop": r["fit_stop"],
            "chisq": r["chisq"],
            "coefficient": r["coefficient"],
            "coefficient_err": r["coefficient_err"],
            "exponent": r["exponent"],
            "exponent_err": r["exponent_err"],
            "constant": constant,
            "constant_err": constant_err,
            "fixedab_constant_err": r["fixedab_constant_err"],
            "best_fit_method": subresults[-1]["fit_method"],
            "best_fit_stop": subresults[-1]["fit_stop"],
            "best_chisq": subresults[-1]["chisq"],
            "best_coefficient": subresults[-1]["coefficient"],
            "best_coefficient_err": subresults[-1]["coefficient_err"],
            "best_exponent": subresults[-1]["exponent"],
            "best_exponent_err": subresults[-1]["exponent_err"],
            "best_constant": best_constant,
            "best_constant_err": best_constant_err,
            "best_fixedab_constant_err": subresults[-1]["fixedab_constant_err"],
            "rel_discrep": constant / best_constant - 1.0,
            "rel_discrep_err": (
                (constant_err / best_constant) ** 2 +
                (constant * best_constant_err / best_constant ** 2) ** 2
            ) ** 0.5,
            "rel_dist": dist / best_constant,
            "rel_dist_err": (
                (r["dist_err"] / best_constant) ** 2 +
                (dist * best_constant_err / best_constant ** 2) ** 2
            ) ** 0.5,
        })
    return results

def gather_fit_data(fit_count, maxfev):
    fn_format = "fits.fit_count={fit_count}_maxfev={maxfev}.txt"
    try:
        return utils.load_table(fn_format.format(**locals()))
    except OSError:
        pass

    d = utils.filter_preferred_ml(utils.load_all())
    d = d[~d["method"].isin(["imsrg[f]+eom[n]"])]
    results = []
    with multiprocessing.Pool(4) as p:
        results = p.map(
            functools.partial(gather_fit_data_inner,
                              fit_count=fit_count,
                              maxfev=maxfev),
            tuple(d.groupby(["label", "interaction", "num_filled",
                             "freq", "method"])))
    d = pd.DataFrame.from_records(itertools.chain(*results))
    print("{} fits failed, out of {}"
          .format((d["fit_method"] == "fixedab").sum(), len(d)))
    # fit_count=5:
    #  maxfev=default: 198 fits failed, out of 2247
    #  maxfev=10k: 40 fits failed, out of 2248
    #  maxfev=100k: 0 fits failed

    cols = """
    interaction
    label
    freq
    num_filled
    method
    best_chisq
    best_coefficient
    best_coefficient_err
    best_constant
    best_constant_err
    best_fixedab_constant_err
    best_exponent
    best_exponent_err
    best_fit_method
    best_fit_stop
    chisq
    coefficient
    coefficient_err
    constant
    constant_err
    fixedab_constant_err
    exponent
    exponent_err
    fit_method
    fit_stop
    rel_discrep
    rel_discrep_err
    rel_dist
    rel_dist_err
    """.split()
    assert len(d.columns) == len(cols)
    utils.save_table(fn_format.format(**locals()), d[cols])
    return d

def parse_bin(s):
    # why Pandas converts bins to strings, no-one knows
    return (tuple(map(float, s[1:-1].split(","))))

def parse_binned_groups(stream):
    return ((parse_bin(k), v) for k, v in stream)

def plot(fit_count=5, log=False, maxfev=0, plot_type="scatter",
         stat="err", hf=False):
    dorig = utils.filter_preferred_ml(utils.load_all())

    d_good = utils.load_table("fits_good.txt")
    d_good = d_good.groupby(
        ["interaction", "label", "freq", "num_filled", "method"]
    ).first()

    d = gather_fit_data(fit_count=fit_count, maxfev=maxfev)

    d = d[d["interaction"] == "normal"]
    # d = d[d["label"] == "add"]
    # d = d[d["method"] == "imsrg"]
    # d = d[d["num_filled"] == 5]
    # d = d[d["freq"] == 0.1]
    d = d[d["fit_method"] != "fixedab"]

    doriggrps = dorig.groupby(["interaction", "label", "freq",
                               "num_filled", "method"])

    d["rel_constant_err"] = d["constant_err"] / d["constant"]
    d["rel_best_constant_err"] = d["best_constant_err"] / d["best_constant"]
    d["label_is_ground"] = d["label"] == "ground"
    d["good"] = d.apply(lambda r: d_good.loc[
        (r["interaction"], r["label"], r["freq"],
         r["num_filled"], r["method"])]["good"], axis=1)
    d["rel_chi"] = d["chisq"]**.5 / d["constant"]
    d["rel_reduced_chi"] = d["rel_chi"] / (fit_count - 3)
    d["rel_best_chisq"] = d["best_chisq"]**.5 / d["best_constant"]
    d["rel_best_reduced_chisq"] = d["rel_best_chisq"] / (fit_count - 3)
    d = d[(d["rel_best_reduced_chisq"] < 1e-6)]
    d["fixedab_with_hf"] = (
        ((d["fit_method"] == "fixedab") ==
         (d["method"].isin(["hf", "hf+qdpt3"]))) |
        (d["fit_method"] == "full")
    )

    color_col = "method"
    bin_transform = TRANSFORM_ID
    bin_transform = TRANSFORM_LOG_ABS

    if color_col in ["exponent", "rel_dist", "rel_constant_err",
                     "rel_best_constant_err", "rel_best_chisq",
                     "rel_chi", "rel_reduced_chi", "chi_ratio"]:
        num_bins = 16
        d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=[color_col])
        binf = bin_transform[0]
        bininvf = bin_transform[1]
        color_bins = pd.cut(binf(abs(d[color_col])), num_bins)
        d["color_bin_start"] = color_bins.map(
            lambda bin: bininvf(parse_bin(bin)[0]))
        d["color_bin_stop"] = color_bins.map(
            lambda bin: bininvf(parse_bin(bin)[1]))
        color_bin_cols = ["color_bin_start", "color_bin_stop"]
    else:
        color_bin_cols = [color_col]
    max_bins = len(d[color_bin_cols[0]].unique())

    fig, ax = plt.subplots()

    def on_pick_event(event):
        x = list(event.artist.get_xdata()[event.ind])[0]
        y = list(event.artist.get_ydata()[event.ind])[0]
        sel = d[(abs(d["x"] - x) <= 1e-20) &
                (abs(d["y"] - y) <= 1e-20)]
        print(sel.transpose().to_csv())
        if len(sel) != 1:
            print('>> not found <<')
            return

        sel = sel.iloc[0]
        grp = doriggrps.get_group((sel["interaction"], sel["label"],
                                   sel["freq"], sel["num_filled"],
                                   sel["method"]))
        fig, ax = plt.subplots(2)
        ax[0].plot(grp["num_shells"], grp["energy"], "x")
        fit_start = sel["fit_stop"] + 1 - fit_count
        ax[0].axvspan(fit_start, sel["fit_stop"], color="#cccccc")
        xs = np.linspace(grp["num_shells"].min(), grp["num_shells"].max())
        ax[0].plot(xs,
                   sel["coefficient"] * xs ** sel["exponent"]
                   + sel["constant"])
        subgrp = grp[grp["num_shells"].between(
            fit_start-0.1, sel["fit_stop"]+0.1)]
        last_constant = sel["constant"]
        last_constant_err = sel["constant_err"]

        def random_weight(count):
            weights = np.zeros(count)
            for i in range(count):
                weights[np.random.randint(0, count)] += 1
            return weights

        p0 = [sel["coefficient"], sel["exponent"], sel["constant"]]
        p = p0
        x = subgrp["num_shells"]
        y = subgrp["energy"]
        constants = []
        constants.append(p[2])

        print(f"x = np.array({list(x)})")
        print(f"y = np.array({list(y)})")

        ax[1].plot(x, (p[0] * x ** p[1] + p[2] - y), "-x")
        ax[1].axhline(0.0, linestyle=":")

        for i in range(10):
            count = len(x)
            weights = random_weight(count) + 1e-99
            if sum(weights > 0.1) <= 3: # can't fit with this few points
                continue
            try:
                p, cov = scipy.optimize.curve_fit(
                    lambda x, a, b, c: a * x ** b + c,
                    x, y,
                    sigma=1.0 / weights ** 0.5,
                    p0=p0, maxfev=100000)
            except RuntimeError as e:
                print(e)
                continue
            chisq = np.average((p[0] * x ** p[1] + p[2] - y) ** 2,
                               weights=weights) * len(x)
            constant = p[2]
            constant_err = cov[2, 2] ** 0.5
            constants.append(p[2])
            last_constant = constant
            last_constant_err = constant_err
        print("result", np.mean(constants), np.std(constants))
        print("rel", np.std(constants) / np.mean(constants))
        ax[0].set_ylim([max(ax[0].get_ylim()[0], 0.0),
                        min(ax[0].get_ylim()[1], np.max(y))])
        ax[0].plot(xs, p[0] * xs ** p[1] + p[2], ":", color="lime")

    fig.canvas.mpl_connect("pick_event", on_pick_event)

    d["x_hessian"] = np.log(d["constant_err"]/d["chisq"]**0.5)

    # hf has unique behaviors (what about mp2?)
    if hf:
        d = d[d["method"] == "hf"]
    else:
        d = d[d["method"] != "hf"]
        if stat == "err":
            d = d[d["x_hessian"] > -0.5]

    d["x"] = d["rel_constant_err"]
    d["y"] = d["rel_discrep"] / d["rel_constant_err"]
    d["x"] = np.log10(d["x"])

    ax.plot(d["x_hessian"], d["constant"], "x")
    return

    if stat == "hessian":
        d["x"] = d["x_hessian"]

    d["y_err"] = d["rel_discrep_err"] / d["rel_discrep"] * d["y"]

    d = d[(d["rel_constant_err"] > 0) & (d["rel_constant_err"] < np.inf)]
    if plot_type == "contour":
        cfg = {
            "err": {
                "nx": 20,
                "ny": 40,
                # ranged = mesh; lim = view
                "lims": {
                    False: { # !hf
                        "xrange": (-7, -1),
                        "yrange": (-50, 50),
                        "xlim": (-6, -2),
                        "ylim": (-40, 40),
                    },
                    True: { # hf
                        "xrange": (-7, -1),
                        "yrange": (-5, 5),
                        "xlim": (-6, -2),
                        "ylim": (-4, 4),
                    },
                },
                "title": ("“actual” discrepancy vs fit uncertainty "
                          "(filtered: Q > -0.5)"),
                "xlabel": r"$\log\left(\frac{\sigma_c}{c}\right)$",
                "ylabel": r"$\frac{\varepsilon}{\sigma_c}$",
            },
            "hessian": {
                "nx": 20,
                "ny": 40,
                # ranged = mesh; lim = view
                "lims": {
                    False: { # !hf
                        "xrange": (-2.0, 8.0),
                        "yrange": (-200, 200),
                        "xlim": (-1.0, 4.0),
                        "ylim": (-150, 150),
                    },
                    True: { # hf
                        "xrange": (-1.0, 2.0),
                        "yrange": (-10, 10),
                        "xlim": (-1.0, 2.0),
                        "ylim": (-10, 10),
                    },
                },
                "title": "spread of “actual” discrepancy",
                "xlabel": (r"$Q = \log\left(\frac{\sigma_c}"
                           r"{\sqrt{\mathtt{RSS}}}\right)$"),
                "ylabel": r"$\frac{\varepsilon}{\sigma_c}$",
            },
        }[stat]
        nx = cfg["nx"]
        ny = cfg["ny"]
        lims = cfg["lims"][hf]
        ax.plot(d["x"], d["y"], "o", markersize=1, picker=3,
                color="white", markeredgewidth=0)
        h, x, y = np.histogram2d(d["x"], d["y"], bins=(nx, ny),
                                 range=(lims["xrange"], lims["yrange"]))
        ch = np.cumsum(h, axis=1)
        z = ch / ch[..., -1, np.newaxis]
        x, y = np.meshgrid(x[:-1], y[1:], indexing="ij")
        levels = np.linspace(-2.0, 2.0, 5)
        levels = 0.5 + 0.5 * scipy.special.erf(2.0 ** -0.5 * levels)
        ax.axhline(-1.0, linestyle="--", color="white", linewidth=0.5)
        ax.axhline(1.0, linestyle="--", color="white", linewidth=0.5)
        ax.contour(x, y, z,
                   levels=levels,
                   colors="white",
                   alpha=0.3)
        cs = ax.contourf(x, y, z,
                         levels=np.linspace(0.0, 1.0, 100),
                         # note: translucent cmaps tend to cause artifacts
                         cmap=CMAP_FOLDED_VIRIDIS,
                         linestyle=":")
        fig.colorbar(cs)
        ax.set_xlabel(cfg["xlabel"])
        ax.set_ylabel(cfg["ylabel"])
        ax.set_xlim(lims["xlim"])
        ax.set_ylim(lims["ylim"])

    elif plot_type == "scatter":

        cmap = utils.CMAP_RAINBOW
        for i, (bin, g) in enumerate(sorted(d.groupby(color_bin_cols))):
            color = cmap(float(i) / max_bins)
            ax.plot(g["x"], g["y"], ".",
                    label=str(bin),
                    color=color,
                    markersize=10,
                    markeredgewidth=0,
                    alpha=0.5,
                    picker=2)
        ax.legend()

    m = min(ax.get_xlim()[1], -ax.get_ylim()[0])

utils.plot_main(__file__, plot, plot)
