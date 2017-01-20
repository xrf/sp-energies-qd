#!/usr/bin/env python3
import os, re, sys
sys.path.insert(1, os.path.join(os.path.dirname(__file__), ".."))
import utils

fn = re.match(r"(.*)-postprocess\.py", __file__).group(1) + ".txt"
d = utils.load_table(fn)
d = utils.check_fun_dep(d,
                        ["interaction", "num_shells", "num_filled", "freq",
                         "method"],
                        {"energy": 2e-5},
                        combiner=utils.rightmost_combiner)
d = d.sort_values(["interaction", "num_shells", "num_filled", "freq",
                   "method", "energy"])
with open(fn, "w") as f:
    f.write("""
# Functional dependencies:
#
#   * (interaction, num_shells, num_filled, freq, method) -> energy
#
"""[1:])
    utils.save_table(f, d)
