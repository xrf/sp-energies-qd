#!/usr/bin/env python3
# Generates ../FigureFiles/fig-compare-*.svg from the data files.
import sys
import matplotlib.pyplot as plt
import pandas as pd
import utils

def plot_compare(suffix, f):
    d = pd.concat(utils.get_ar_energies_for_v(suffix))
    sys.stderr.write("Saving ...\n".format(**locals()))
    sys.stderr.flush()
    for freq in [0.28, 1.0]:
        for num_filled in [2, 3, 4]:
            num_particles = num_filled * (num_filled + 1)
            for label in ["add", "rm"]:
                fig, axs = plt.subplots(2, 1)
                fig.set_size_inches(8, 10)
                ml = (num_filled + (label != "add")) % 2
                for zoomed in [False, True]:
                    ax = axs[int(zoomed)]
                    # zoom in to the part where it's nearly converged
                    num_shells_start = {
                        2: 7,
                        3: 9,
                        4: 11,
                    }[num_filled] if zoomed else 5
                    for method in d["method"].unique():
                        xs = []
                        ys = []
                        for num_shells in range(num_shells_start, 16):
                            case = d[(d["num_shells"] == num_shells) &
                                     (d["num_filled"] == num_filled) &
                                     (d["freq"] == freq) &
                                     (d["ml"] == ml) &
                                     (d["label"] == label)]
                            try:
                                energy, = case[case["method"] == method].energy
                            except ValueError:
                                continue
                            xs.append(num_shells)
                            ys.append(energy)
                        kwargs = {
                            "qdpt": {
                                "color": "#39b237",
                            },
                            "eom": {
                                "color": "#841a0b",
                            },
                            "eom_f": {
                                "color": "#c6cc26",
                                "linestyle": "--",
                            },
                            "eom_quads": {
                                "color": "#a825bc",
                            },
                            "cc": {
                                "color": "#1351c4",
                            },
                        }[method]
                        ax.plot(xs, ys, "-x", label=method,
                                linewidth=1.5, **kwargs)
                        if method == "qdpt":
                            slope = ((ys[-1] - ys[-2]) / (xs[-1] - xs[-2]) / ys[-1])
                    if not zoomed:
                        f.write(
                            "{num_particles}\t{freq}\t{label}\t{suffix!r}\t{slope}\n"
                            .format(**locals()))
                        ax.set_title("{} energy for {} particles (ML = {}, omega = {})\n"
                                     "{} (rel_slope = {:.6f})"
                                     .format({"add": "addition",
                                              "rm": "removal"}[label],
                                             num_particles, ml, freq, suffix, slope), y=1.08)
                        ax.legend()
                    ax.set_xlabel("number_of_shells")
                    ax.set_ylabel("energy")
                fn = ("../FigureFiles/fig-compare-"
                      "{num_particles}-{freq}-{label}-{ml}{suffix}.svg"
                      .format(**locals()))
                sys.stderr.write(fn + "\n")
                sys.stderr.flush()
                fig.savefig(fn)
                plt.close(fig)

utils.init(__file__)

with open("compare_rel_slopes.txt", "w") as f:
    f.write("# num_particles\tfreq\tlabel\tinteraction\tslope\n")
    for suffix in ["", "_sigmaA=0.5_sigmaB=4.0", "_sigmaA=0.5", "_sigmaB=4.0"]:
        plot_compare(suffix, f)
