#!/usr/bin/env python3
import functools, random
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize
import scipy.special
import fits, utils

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

def parse_bin(s):
    # why Pandas converts bins to strings, no-one knows
    return (tuple(map(float, s[1:-1].split(","))))

def parse_binned_groups(stream):
    return ((parse_bin(k), v) for k, v in stream)

def is_good(d_good, r):
    key = r["interaction"], r["label"], r["freq"], r["num_filled"], r["method"]
    try:
        return d_good.loc[key]["good"]
    except KeyError:
        # unknown: default to False
        return False

def plot(fit_count=fits.DEFAULT_FIT_COUNT, log=False,
         maxfev=fits.DEFAULT_MAXFEV, plot_type="scatter",
         stat="err", hf=False):
    dorig = utils.filter_preferred_ml(utils.load_all())

    d_good = utils.load_table("fits_good.txt")
    d_good = d_good.groupby(
        ["interaction", "label", "freq", "num_filled", "method"]
    ).first()

    d = fits.load_predictive_data(fit_count=fit_count, maxfev=maxfev)

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
    d["good"] = d.apply(functools.partial(is_good, d_good), axis=1)
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

    d["quality"] = np.log10(d["constant_err"]/d["chisq"]**0.5)

    # hf has unique behaviors (what about mp2?)
    if hf:
        d = d[d["method"] == "hf"]
    else:
        d = d[d["method"] != "hf"]
        if stat == "err":
            d = d[d["quality"] > 0]

    d["y"] = d["rel_discrep"] / d["rel_constant_err"]
    if stat == "hessian":
        d["x"] = d["quality"]
    elif stat == "err":
        d["x"] = np.log10(d["rel_constant_err"])
    else:
        assert False

    d["y_err"] = d["rel_discrep_err"] / d["rel_discrep"] * d["y"]

    d = d[(d["rel_constant_err"] > 0) & (d["rel_constant_err"] < np.inf)]
    if plot_type == "contour":

        if hf:
            hf_suffix = "HF"
        else:
            hf_suffix = "non-HF"
        # ranged = mesh; lim = view
        if stat == "err":
            if hf:
                nx = 20
                ny = 20
                xrange = (-7, -1)
                yrange = (-5, 5)
                xlim = (-6, -2)
                ylim = (-4, 4)
                title = ("discrepancy vs fit uncertainty "
                         f"({hf_suffix})")
            else:
                nx = 20
                ny = 40
                xrange = (-7, -1)
                yrange = (-50, 50)
                xlim = (-6, -2)
                ylim = (-40, 40)
                title = ("discrepancy vs fit uncertainty "
                         f"(filtered: Q > 0, {hf_suffix})")
            xlabel = r"$\log_{10}\left(\frac{\sigma_c}{c}\right)$"
            ylabel = r"$\frac{\varepsilon}{\sigma_c}$"
        elif stat == "hessian":
            if hf:
                nx = 20
                ny = 40
                xrange = (-0.6, 1.5)
                yrange = (-10, 10)
                xlim = (-0.4, 1.1)
                ylim = (-10, 10)
            else:
                nx = 20
                ny = 40
                xrange = (-0.6, 2.1)
                yrange = (-200, 200)
                xlim = (-0.4, 1.7)
                ylim = (-150, 150)
            title = ("spread of “actual” discrepancy vs quality "
                     f"({hf_suffix})")
            xlabel = (r"$Q = \log_{10}\left(\frac{\sigma_c}"
                       r"{\sqrt{\mathtt{RSS}}}\right)$")
            ylabel = r"$\frac{\varepsilon}{\sigma_c}$"
        else:
            assert False
        ax.plot(d["x"], d["y"], "o", markersize=1, picker=3,
                color="white", markeredgewidth=0)
        dx = (xrange[1] - xrange[0]) / (nx - 1)
        h, x, y = np.histogram2d(d["x"], d["y"], bins=(nx, ny - 1),
                                 range=((xrange[0] - 0.5 * dx,
                                         xrange[1] + 0.5 * dx),
                                        yrange))
        ch = np.concatenate([np.zeros((nx, 1)),
                             np.cumsum(h, axis=1)],
                            axis=1)
        z = ch / ch[..., -1, np.newaxis]
        x, y = np.meshgrid(0.5 * (x[1:] + x[:-1]), y, indexing="ij")
        levels = np.linspace(-2.0, 2.0, 5)
        levels = 0.5 + 0.5 * scipy.special.erf(2.0 ** -0.5 * levels)
        ax.axhline(-1.0, linestyle="--", color="white", linewidth=0.5)
        ax.axhline(1.0, linestyle="--", color="white", linewidth=0.5)
        ax.contour(x, y, z,
                   levels=levels,
                   colors="white",
                   alpha=0.3)
        cs = ax.contourf(x, y, z,
                         levels=np.linspace(0.0, 1.0, 300),
                         # note: translucent cmaps tend to cause artifacts
                         cmap=CMAP_FOLDED_VIRIDIS,
                         linestyle=":")
        for c in cs.collections: # http://stackoverflow.com/a/32911283
            c.set_edgecolor("face")
        fig.colorbar(cs)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

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

    fn = f"fit-predictiveness-{plot_type}-{fit_count}-{stat}-{hf}"
    fig.tight_layout()
    utils.savefig(fig, fn)

def main():
    maxfev = fits.DEFAULT_MAXFEV
    fit_count = fits.DEFAULT_FIT_COUNT
    for stat in ["err", "hessian"]:
        for hf in [False, True]:
            plot(fit_count=fit_count,
                 maxfev=maxfev,
                 plot_type="contour",
                 stat=stat,
                 hf=hf)

utils.plot_main(__file__, plot, main)
