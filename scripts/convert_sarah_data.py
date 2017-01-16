#!/usr/bin/env python3
import argparse, re, string
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("in_file")
p.add_argument("out_file")
kwargs = vars(p.parse_args())

pos_to_col = {
    0: "num_shells",
    1: "num_shells",
    8: "energy_hf",
    24: "energy_imsrg",
    25: "energy_imsrg",
    40: "energy_ccsd",
    56: "energy_fci",
}

col_types = {
    "num_shells": int,
    "energy_hf": float,
    "energy_imsrg": float,
    "energy_ccsd": float,
    "energy_fci": float,
}

data = []
with open(kwargs["in_file"]) as f:
    for line in f:
        if not line.strip():
            continue

        line = line.expandtabs()

        m = re.match(r"##### (\d+) ([\d.]+) #####\s*$", line)
        if m:
            num_particles, freq = m.groups()
            continue

        if re.match("no convergence", line):
            continue

        if not re.match("(\s*([\d.]+|nc))*\s*$", line):
            raise Exception("unrecognized line: {!r}".format(line))

        rec = {}
        for m in re.finditer("\s*([\d.]+|nc)", line):
            value, = m.groups()
            if value == "nc":
                value = "nan"
            col = pos_to_col[m.start(1)]
            rec[col] = col_types[col](value)

        for k, v in rec.items():
            if k == "num_shells":
                continue
            method = re.match("energy_(.*)", k).group(1)
            data.append({
                "num_particles": num_particles,
                "freq": freq,
                "num_shells": rec["num_shells"],
                "method": method,
                "energy": v,
            })

pd.DataFrame.from_dict(data)[[
    "num_particles",
    "freq",
    "num_shells",
    "method",
    "energy",
]].to_csv(kwargs["out_file"], sep=" ", na_rep="nan", index=False)
