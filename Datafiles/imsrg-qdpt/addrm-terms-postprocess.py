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
                         "method", "p", "term_id"],
                        {"correction": 1e-7},
                        combiner=utils.rightmost_combiner)
d = d.sort_values(["interaction", "num_shells", "num_filled", "freq",
                   "method", "p", "term_id", "correction"])
with open(fn, "w") as f:
    f.write("""
# term_ids 3 and 4: QDPT2
# term_ids 5 to 22: QDPT3
#
# Functional dependencies:
#
#   * (num_shells, num_filled, freq, method, p, term_id) -> correction
#
"""[1:])
    utils.save_table(f, d)
