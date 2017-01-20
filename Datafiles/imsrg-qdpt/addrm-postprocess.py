#!/usr/bin/env python3
import os, re, sys
sys.path.insert(1, os.path.join(os.path.dirname(__file__), ".."))
import utils

fn = re.match(r"(.*)-postprocess\.py", __file__).group(1) + ".txt"
d = utils.load_table(fn)
# canonicalization can introduce duplicates, in addition to whatever
# duplicates that already exist in the file
d["p"] = d["p"].map(utils.canonicalize_p)
d = utils.check_fun_dep(d,
                        ["interaction", "num_shells", "num_filled", "freq",
                         "method", "p", "iter"],
                        {"energy": 1e-6},
                        combiner=utils.rightmost_combiner)
d = d.sort_values(["interaction", "num_shells", "num_filled", "freq",
                   "method", "p", "iter", "energy"])
with open(fn, "w") as f:
    f.write("""
# iter = 0: without QDPT3 correction
# iter = 1: includes QDPT3 correction
#
# Functional dependencies:
#
#   * (interaction, num_shells, num_filled, freq, method, p) -> correction
#
"""[1:])
    utils.save_table(f, d)
