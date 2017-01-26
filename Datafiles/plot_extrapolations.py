#!/usr/bin/env python3
import json, logging, os
import utils
import pandas as pd
import matplotlib.pyplot as plt

def convert_fit_json_to_dataframe(json_data):
    data = []
    for entry in json_data:
        entry = dict(entry)
        fit_result = entry.pop("fit_result")

        r = fit_result.get("full")
        if r:
            entry["fit_method"] = "full"
            entry["constant"] = r["constant"]
            entry["constant_err"] = r["constant_err"]
            data.append(entry)
            continue

        r = fit_result.get("fixedab")
        if "fixedab" in fit_result:
            entry["fit_method"] = "fixedab"
            entry["constant"] = r["constant"]
            # TMP: bugfix:
            if isinstance(r["constant_err"], list):
                r["constant_err"] = r["constant_err"][0]
            entry["constant_err"] = r["constant_err"]
            data.append(entry)
            continue

        logging.warning(f"no fit result available for: {entry}")
    return pd.DataFrame.from_records(data)

os.chdir(os.path.dirname(__file__))
with open("fit_results.json") as f:
    j = json.load(f)
d = convert_fit_json_to_dataframe(j)

label = "ground"
interaction = "normal"
if True:
    num_filled = 3
    gs = d.groupby(["label", "interaction", "num_filled"])
    d = gs.get_group((label, interaction, num_filled))
    fig, ax = plt.subplots()
    for method, g in d.groupby("method"):
        x = g["freq"]
        y = g["constant"]
        y_err = g["constant_err"]
        ax.plot(x, y, label=method, marker="x")
        ax.fill_between(x, y - y_err, y + y_err, label="", alpha=0.3)
else:
    freq = 1.0
    gs = d.groupby(["label", "interaction", "freq"])
    d = gs.get_group((label, interaction, freq))
    fig, ax = plt.subplots()
    for method, g in d.groupby("method"):
        x = g["num_filled"]
        y = g["constant"] / g["num_filled"]
        y_err = g["constant_err"] / g["num_filled"]
        ax.plot(x, y, label=method, marker="x")
        ax.fill_between(x, y - y_err, y + y_err, label="", alpha=0.3)
ax.legend()
plt.show()
