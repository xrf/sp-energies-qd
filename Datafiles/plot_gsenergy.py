#!/usr/bin/env python3
import argparse, json, os, scipy, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utils

FIT_TYPES = ["loglog", "semilog"]

def plot(interactive):
    if interactive:
        plt.ion()

    d = pd.read_csv("imsrg-qdpt/dat_gsenergy.txt",
                    header=0, index_col=False,
                    delim_whitespace=True)
    d = d[d["num_filled"] == 3]
    d = d[d["freq"] == 1.0]
    d = d[d["freq"] == 1.0]
    fit_type = "loglog"
    fit_points = 7

    for [num_filled, freq], gg in d.groupby(["num_filled", "freq"]):
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 10)
        for method, g in gg.groupby(["method"]):
            g = g.sort_values(["num_shells"])
            g["d_num_shells"] = g["num_shells"].diff()
            g["d_energy"] = g["energy"].diff()
            g["d_energy_d_num_shells"] = g["energy"].diff()
            #g = g.dropna()

            def guess_constant(a, b):
                return g["energy"].iat[-1] - a * g["num_shells"].iat[-1] ** b

            x = g["num_shells"] - g["d_num_shells"] / 2.0
            y = abs(g["d_energy"] / g["d_num_shells"])
            fit = utils.fit_change(fit_type=fit_type,
                                   fit_points=fit_points,
                                   x=x,
                                   y=y,
                                   num_filled=num_filled,
                                   freq=freq,
                                   method=method)
            sys.stdout.write("{}\n\n".format(
                json.dumps(utils.sanitize_json(fit["extra"]),
                           **utils.JSON_PRETTY)))
            sys.stdout.flush()
            guess_c = guess_constant(fit["extra"]["coefficient"],
                                     fit["extra"]["exponent"])
            fit_label = ("{method} deriv fit ("
                         "E_gnd = ({coefficient:.3g}±{coefficient_err:.2g}) * K ** "
                         "({exponent:.3g}±{exponent_err:.2g}) + "
                         "(≈{constant:.5g})"
                         .format(constant=guess_c, **fit["extra"]))

            def f(x, a, b, c):
                return a * x ** b + c
            # parameter guesses using the loglog fit (NOT the nonlinear power
            # fit, because the nonlinear fit sometimes gets "stuck" when the
            # data is bad)
            b = fit["guess"][0] + 1.0
            a = np.exp(fit["guess"][1]) / b
            c = guess_constant(a, b)
            p, var_p = scipy.optimize.curve_fit(
                f,
                g["num_shells"][-(fit_points + 1):],
                g["energy"][-(fit_points + 1):],
                p0=(a, b, c),
            )
            fit2_results = {
                "num_filled": num_filled,
                "freq": freq,
                "method": method,
                "coefficient": p[0],
                "coefficient_err": var_p[0, 0] ** 0.5,
                "exponent": p[1],
                "exponent_err": var_p[1, 1] ** 0.5,
                "constant": p[2],
                "constant_err": var_p[2, 2] ** 0.5,
            }
            sys.stdout.write("{}\n\n".format(
                json.dumps(utils.sanitize_json(fit2_results),
                           **utils.JSON_PRETTY)))
            sys.stdout.flush()
            fit2_label = (
                "{method} deriv fit ("
                "E_gnd = ({coefficient:.3g}±{coefficient_err:.2g}) * K ** "
                "({exponent:.3g}±{exponent_err:.2g}) + "
                "({constant:.5g}±{constant_err:.2g})"
                .format(**fit2_results))
            ax.loglog(x, y, "o", label=method,
                      color=utils.GS_METHOD_COLOR[method])
            ax.loglog(fit["x"], fit["y"], "-", label=fit_label,
                      color=utils.GS_METHOD_COLOR[method])
            ax.loglog(x, abs(p[0] * p[1]) * x ** (p[1] - 1), ":",
                      label=fit2_label, color=utils.GS_METHOD_COLOR[method])
        ax.legend()
        ax.set_xlabel("K (number of shells)")
        ax.set_ylabel("|ΔE_gnd/ΔK|")

    if not interactive:
        fn = "fig-gsenergy2-change-{num_filled}-{freq}.svg".format(**locals())
        fig.savefig(fn)
        plt.close(fig)
        sys.stderr.write("// Figure saved to: {}\n\n".format(fn))
        sys.stderr.flush()

    for [num_filled, freq], gg in d.groupby(["num_filled", "freq"]):
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 10)
        for method, g in gg.groupby(["method"]):
            ax.plot(
                g["num_shells"],
                g["energy"] / (num_filled * (num_filled + 1)),
                label=method,
            )
        ax.legend()
        ax.set_xlabel("K (number of shells)")
        ax.set_ylabel("E_gnd/N (energy per particle)")

    if not interactive:
        fn = "fig-gsenergy2-{num_filled}-{freq}.svg".format(**locals())
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
    plot(**kwargs)

main()
