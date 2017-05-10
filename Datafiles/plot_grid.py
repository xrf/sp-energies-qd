#!/usr/bin/env python3
import itertools, os
import matplotlib.lines
import matplotlib.gridspec
import matplotlib.pyplot as plt
import utils

def main(label):
    num_shells_range = {
        6: [3.5, 14.5],
        12: [4.5, 16.5],
        20: [6.5, 16.5],
        30: [9.5, 16.5],
        42: [11.5, 20.5],
        56: [14.5, 20.5],
    }
    num_particles_ticks = [6, 30, 56]
    LEGEND_FACET_IDX = 2

    if label == "ground":
        freq_ticks = [0.28, 1.0]
        methods = ["mp2", "imsrg", "ccsd", "fci"]
        fn = "../Manuscript/fig-gs2.pdf"
    else:
        freq_ticks = [0.28, 1.0]        # 0.28 needed to show DMC results
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

    has_dmc = False
    width = 6.5
    height = 3
    base_markersize = 5.0
    fig = plt.figure(figsize=(width, height))
    gs = matplotlib.gridspec.GridSpec(len(facet_y["ticks"]),
                                      len(facet_x["ticks"]))
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
            has_dmc = True
            ax.axhline(sel_dmc["energy"],
                       color=utils.METHOD_COLOR["dmc"],
                       linestyle=utils.DMC_LINESTYLE)

        for method, d in d.groupby("method"):
            if method not in methods:
                continue

            d = d[(d["num_shells"] >= num_shells_range[num_particles][0]) &
                  (d["num_shells"] <= num_shells_range[num_particles][1])]
            d = d.sort_values(["num_shells"])
            marker = utils.METHOD_MARKER[method]
            ax.plot(d["num_shells"], d["energy"],
                    marker=marker,
                    markerfacecolor="none",
                    markersize=(utils.MARKERSIZE_CORRECTION.get(marker, 1.0) *
                                base_markersize),
                    color=utils.METHOD_COLOR[method])

        title = ("${} = {}$\n${} = {}$"
                 .format(facet_x["symbol"], fct_x,
                         facet_y["symbol"], fct_y))
        ax.text(1.0 - 2.65 / width, 0.6, title,
                color="#505050",
                horizontalalignment="left",
                transform=ax.transAxes,
                fontsize=12)
        ax.set_xlim(num_shells_range[num_particles])
    fig.text(0.5, 0.05 / height, "$K$ (number of shells)",
             horizontalalignment="center",
             verticalalignment="bottom",
             transform=ax.transAxes)
    fig.text(0.05 / width, 0.5, "${}$".format(utils.ENERGY_SYMBOL[label]),
             horizontalalignment="left",
             verticalalignment="center",
             transform=ax.transAxes,
             rotation="vertical")

    # phantom lines to configure the legend
    markersize = (utils.MARKERSIZE_CORRECTION.get(marker, 1.0) *
                  base_markersize)
    lines = [matplotlib.lines.Line2D([], [],
                                     marker=utils.METHOD_MARKER[method],
                                     markersize=markersize,
                                     color=utils.METHOD_COLOR[method],
                                     label=utils.METHOD_LABEL[method])
             for method in methods]
    if has_dmc:
        lines.append(matplotlib.lines.Line2D([], [],
                                             linestyle=utils.DMC_LINESTYLE,
                                             color=utils.METHOD_COLOR["dmc"],
                                             label=utils.METHOD_LABEL["dmc"]))
    ax = fig.add_subplot(gs[LEGEND_FACET_IDX])
    ax.axis("off")
    ax.legend(handles=lines, loc="center", frameon=False,
              bbox_to_anchor=(0.5, 0.5 - 0.1 / height))

    left_margin = 0.1
    right_margin = 0.0
    top_margin = 0.0
    bottom_margin = 0.1
    gs.tight_layout(fig, rect=[left_margin / width,
                               bottom_margin / height,
                               1.0 - right_margin / width,
                               1.0 - top_margin / height])
    utils.savefig(fig, fn=fn)

with utils.plot(__file__):
    main("ground")
    main("add")
    main("rm")
