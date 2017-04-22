import utils

DEFAULT_FIT_COUNT = 5
DEFAULT_MAXFEV = 100000

class NoMatchButBestExistsError(KeyError):
    pass

def load_fit_data(label, fit_count=DEFAULT_FIT_COUNT, maxfev=DEFAULT_MAXFEV):
    d = utils.load_table(f"fits.fit_count={fit_count}_maxfev={maxfev}.txt")
    d = d[d["interaction"] == "normal"]
    d = d[d["label"] == label]
    d["num_particles"] = d["num_filled"] * (d["num_filled"] + 1)
    d_best = d.groupby(["num_particles", "freq", "method"]).first()
    d = d.set_index(["fit_stop", "num_particles", "freq", "method"])
    return {"fit": d, "fit_best": d_best}

def get_fit_params(fit_data, num_shells, num_particles, freq, method):
    '''Raises KeyError if the data cannot be found.'''
    ks = ["coefficient",
          "coefficient_err",
          "exponent",
          "exponent_err",
          "constant",
          "constant_err"]
    try:
        r = fit_data["fit"].loc[(num_shells, num_particles, freq, method)]
    except KeyError:
        r = None
    else:
        return dict((k, r[k]) for k in ks)
    if r is None:
        # some of the fit_stops are "best"
        # so they aren't included directly
        try:
            r = fit_data["fit_best"].loc[(num_particles, freq, method)]
        except KeyError:
            raise KeyError(r) from None
        if r["best_fit_stop"] != num_shells:
            raise NoMatchButBestExistsError((num_particles, freq, method),
                                            {"best": r["best_fit_stop"]})
        return dict((k, r["best_" + k]) for k in ks)
