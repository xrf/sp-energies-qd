#!/usr/bin/env python3
import itertools, os
import matplotlib.lines
import matplotlib.gridspec
import matplotlib.pyplot as plt
import utils

def main(label):
    num_shells_range = {
        6: [3.5, 15.5],
        12: [4.5, 15.5],
        20: [6.5, 15.5],
        30: [9.5, 20.5],
        42: [11.5, 20.5],
        56: [14.5, 20.5],
    }
    num_particles_ticks = [6, 12, 20, 30, 42, 56]
    freq_ticks = [0.1, 0.28, 1.0]

    if label == "ground":
        methods = ["mp2", "imsrg", "ccsd", "fci"]
        fn = "../Manuscript/fig-gs2.pdf"
    else:
        fn = "../Manuscript/fig-{}2.pdf".format(label)
        methods = ["imsrg+qdpt3", "imsrg+eom", "ccsd+eom"]

    d = utils.load_all()
    d_dmc = utils.load_all_dmc()
    d = utils.filter_preferred_ml(d)
    del d["ml"]
    d = d[d["interaction"] == "normal"]
    del d["interaction"]
    d = d[d["label"] == label]
    d_dmc = d_dmc[d_dmc["label"] == label]
    del d["label"]
    del d_dmc["label"]
    del d["num_filled"]
    d_dmc = d_dmc.set_index(["num_particles", "freq"])

    facet_x = {
        "col": "num_particles",
        "symbol": "N",
        "ticks": num_particles_ticks,
    }
    facet_y = {
        "col": "freq",
        "symbol": "\omega",
        "ticks": freq_ticks,
    }
    facet_to_num_particles_freq = lambda x, y: (x, y)

    base_markersize = 3.0
    linewidth = 1.0
    fig = plt.figure(figsize=(9, 4))
    gs = matplotlib.gridspec.GridSpec(len(facet_y["ticks"]),
                                      len(facet_x["ticks"]))
    #d = d.set_index(["num_shells", "num_particles", "freq", "method"])
    LEGEND_FACET_IDX = 5
    for (fct_x, fct_y), d in d.groupby([facet_x["col"], facet_y["col"]]):
        num_particles, freq = facet_to_num_particles_freq(fct_x, fct_y)
        try:
            gi = (facet_x["ticks"].index(fct_x) +
                  facet_y["ticks"].index(fct_y) * len(facet_x["ticks"]))
        except ValueError:
            continue
        if gi == LEGEND_FACET_IDX:
            continue
        ax = fig.add_subplot(gs[gi])

        try:
            sel_dmc = d_dmc.loc[(num_particles, freq)]
        except KeyError:
            pass
        else:
            ax.axhline(sel_dmc["energy"],
                       color=utils.METHOD_COLOR["dmc"],
                       linestyle=utils.DMC_LINESTYLE,
                       linewidth=linewidth)

        for method, d in d.groupby("method"):
            if method not in methods:
                continue

            d = d[(d["num_shells"] >= num_shells_range[num_particles][0]) &
                  (d["num_shells"] <= num_shells_range[num_particles][1])]
            d = d.sort_values(["num_shells"])
            marker = utils.METHOD_MARKER[method]
            ax.plot(d["num_shells"], d["energy"],
                    linewidth=linewidth,
                    marker=marker,
                    markersize=(utils.MARKERSIZE_CORRECTION.get(marker, 1.0) *
                                base_markersize),
                    color=utils.METHOD_COLOR[method])

            title = ("${}={}$\n${}={}$"
                     .format(facet_x["symbol"], fct_x,
                             facet_y["symbol"], fct_y))
            ax.text(.4, .7, title,
                    color="#666666",
                    horizontalalignment="left",
                    transform=ax.transAxes)
        ax.set_xlim(num_shells_range[num_particles])

    # phantom lines to configure the legend
    markersize = (utils.MARKERSIZE_CORRECTION.get(marker, 1.0) *
                  base_markersize)
    lines = ([matplotlib.lines.Line2D([], [],
                                      linewidth=linewidth,
                                      marker=utils.METHOD_MARKER[method],
                                      markersize=markersize,
                                      color=utils.METHOD_COLOR[method],
                                      label=utils.METHOD_LABEL[method])
             for method in methods] +
             [matplotlib.lines.Line2D([], [],
                                      linewidth=linewidth,
                                      linestyle=utils.DMC_LINESTYLE,
                                      color=utils.METHOD_COLOR["dmc"],
                                      label=utils.METHOD_LABEL["dmc"])])
    ax = fig.add_subplot(gs[LEGEND_FACET_IDX])
    ax.axis("off")
    ax.legend(handles=lines)

    fig.tight_layout()
    utils.savefig(fig, fn=fn)

with utils.plot(__file__):
    matplotlib.rcParams["axes.titlesize"] = "medium"
    matplotlib.rcParams["font.size"] = 8
    main("ground")
    main("add")
    main("rm")
