#!/usr/bin/env python
"""Averaged background bubble-rate plot vs. Seitz threshold.

Reads the refactor exposure tables written by run_exposures.py
(``<stem>_exposures_mix.txt`` in ``exposures-refactor``), converts each fitted
nucleation lifetime tau into a rate, and combines runs taken at the same
operating temperature. Points are plotted against the Seitz threshold [keV]
(read straight from the file), one series per temperature.

Only the four SeitzModel columns in the input encode temperature; the temperature
value itself is re-derived from each config's run dates via temperature_K() (the
same date-cut rule used by run_exposures.py).

Usage:
    python rate_plot.py                 # background configs -> plots_07_01_26
    python rate_plot.py --no-plots      # just print the averaged-rate table
    python rate_plot.py --prefix "Cold Cs" --name cold_cs_rate
"""

import argparse
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from configs import CONFIGS, sanitize

# ----------------------------------------------------------------------------- #
# Defaults / constants
# ----------------------------------------------------------------------------- #
DEFAULT_INDIR = "/exp/e961/data/users/gputnam/exposures-refactor"
DEFAULT_PLOTDIR = "/exp/e961/app/users/gputnam/LAr10Ana/plots_07_01_26"
DEFAULT_PREFIX = "Background"     # config-title prefix selecting the data scope
DEFAULT_NAME = "background_rate"

COLORS = ["black", "lightslategray"]
SEC_PER_HOUR = 3600.0

# Quality cuts: drop degenerate low-statistics fits
# (tau rails to a tiny boundary value) and points with implausible fit precision.
DEFAULT_MIN_TAU = 10.0        # s
DEFAULT_MIN_FRACERR = 0.05    # tau_err / tau
DEFAULT_MAX_FRACERR = 0.30

# Input column indices (see OUT_HEADER in run_exposures.py).
C_PRESSURE = 0
C_TAU = 1
C_TAU_ERR = 2
C_SEITZ = 8

def temperature_K(runs):
    """Operating temperature (K) for a configuration from its run dates."""
    dates = [int(r.split("_")[0]) for r in runs]
    before = [d < TEMP_DATE_CUT for d in dates]
    if all(before):
        return TEMP_BEFORE_K
    if not any(before):
        return TEMP_AFTER_K
    # Configurations are date-contiguous and never straddle the cut, but guard
    # anyway: use the value for the earliest run and warn.
    print("  WARNING: config straddles temperature cut %d; using earliest-run value"
          % TEMP_DATE_CUT)
    return TEMP_BEFORE_K if min(dates) < TEMP_DATE_CUT else TEMP_AFTER_K

# ----------------------------------------------------------------------------- #
# Data loading / rate computation
# ----------------------------------------------------------------------------- #
def load_points(titles, indir, min_tau, min_fracerr, max_fracerr):
    """Collect per-run rate points, tagged by temperature.

    Returns a dict {temp_K: [ (pressure, seitz, rate, rate_err), ... ]} where rate
    is in bubbles/hour.  Quality cuts (matching the source notebook) drop failed
    fits and degenerate low-statistics fits: tau > min_tau and
    min_fracerr < tau_err/tau < max_fracerr.  Files without a matching
    _exposures_mix.txt are skipped with a warning.
    """
    points = {}
    for title in titles:
        stem = sanitize(title)
        path = os.path.join(indir, "%s_exposures_mix.txt" % stem)
        if not os.path.exists(path):
            print("  WARNING: no file for %r (%s); skipping" % (title, path))
            continue

        data = np.atleast_2d(np.loadtxt(path))
        if data.size == 0:
            print("  WARNING: empty file for %r; skipping" % title)
            continue

        temp = temperature_K(CONFIGS[title])

        tau = data[:, C_TAU]
        tau_err = data[:, C_TAU_ERR]
        with np.errstate(divide="ignore", invalid="ignore"):
            fracerr = tau_err / tau
        good = ((tau > min_tau) & (tau_err > 0)
                & (fracerr > min_fracerr) & (fracerr < max_fracerr))

        rate = SEC_PER_HOUR / tau[good]
        rate_err = SEC_PER_HOUR * tau_err[good] / tau[good] ** 2
        press = data[good, C_PRESSURE]
        seitz = data[good, C_SEITZ]

        bucket = points.setdefault(temp, [])
        for p, q, r, re in zip(press, seitz, rate, rate_err):
            bucket.append((p, q, r, re))

        print("  %-24s T=%.1f K  %d valid setpoint(s)"
              % (title, temp, int(good.sum())))
    return points


def average_by_setpoint(rows):
    """Inverse-variance weighted mean of per-run rates, grouped by pressure.

    ``rows`` is a list of (pressure, seitz, rate, rate_err).  Returns lists
    (pressure, seitz, rate_avg, rate_avg_err, n_runs) sorted by Seitz threshold.
    """
    groups = {}
    for p, q, r, re in rows:
        groups.setdefault(round(p, 2), []).append((q, r, re))

    out = []
    for p, items in groups.items():
        q = np.array([i[0] for i in items])
        r = np.array([i[1] for i in items])
        re = np.array([i[2] for i in items])
        w = 1.0 / re ** 2
        rate_avg = np.sum(w * r) / np.sum(w)
        rate_avg_err = np.sqrt(1.0 / np.sum(w))
        out.append((p, float(np.mean(q)), rate_avg, rate_avg_err, len(items)))

    out.sort(key=lambda row: row[1])   # by Seitz threshold
    return out


# ----------------------------------------------------------------------------- #
# Output
# ----------------------------------------------------------------------------- #
def print_table(averaged):
    """Print the averaged-rate table to stdout (one block per temperature)."""
    print("\nAveraged background rate (inverse-variance weighted per setpoint):")
    header = "  %-8s %-14s %-12s %-22s %s" % (
        "T[K]", "Pressure[bara]", "Seitz[keV]", "Rate[bubbles/hour]", "N runs")
    for temp in sorted(averaged):
        print(header)
        for p, q, r, re, n in averaged[temp]:
            print("  %-8.1f %-14.2f %-12.3f %-22s %d"
                  % (temp, p, q, "%.3f +/- %.3f" % (r, re), n))
        print("")


def make_plot(averaged, plotdir, name):
    """Errorbar plot of averaged rate vs. Seitz threshold, one series per T."""
    fig = plt.figure()
    for i, temp in enumerate(sorted(averaged)):
        rows = averaged[temp]
        seitz = [row[1] for row in rows]
        rate = [row[2] for row in rows]
        rate_err = [row[3] for row in rows]
        plt.errorbar(seitz, rate, rate_err, linestyle="none", marker=".",
                     label="$%.1f\\,$K" % temp, color=COLORS[i % len(COLORS)])

    l = plt.legend(title="Background Rate\nat Temperature:")
    l.get_title().set_ha("center")
    plt.xlabel("Seitz Threshold [keV]")
    plt.ylabel("Bubble Rate [bubbles / hour]")
    plt.text(0.975, 0.55, "SBC Preliminary", fontsize=16, fontweight="bold",
             horizontalalignment="right", color="tab:purple",
             transform=plt.gca().transAxes)

    os.makedirs(plotdir, exist_ok=True)
    for ext in ("pdf", "png"):
        out = os.path.join(plotdir, "%s.%s" % (name, ext))
        fig.savefig(out)
        print("  wrote %s" % out)
    plt.close(fig)


# ----------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=DEFAULT_INDIR)
    ap.add_argument("--plotdir", default=DEFAULT_PLOTDIR)
    ap.add_argument("--prefix", default=DEFAULT_PREFIX,
                    help="config-title prefix to select (default %r)" % DEFAULT_PREFIX)
    ap.add_argument("--name", default=DEFAULT_NAME,
                    help="output filename stem (default %r)" % DEFAULT_NAME)
    ap.add_argument("--min-tau", type=float, default=DEFAULT_MIN_TAU,
                    help="drop fits with tau <= this many s (default %.1f)"
                         % DEFAULT_MIN_TAU)
    ap.add_argument("--min-fracerr", type=float, default=DEFAULT_MIN_FRACERR,
                    help="drop fits with tau_err/tau <= this (default %.2f)"
                         % DEFAULT_MIN_FRACERR)
    ap.add_argument("--max-fracerr", type=float, default=DEFAULT_MAX_FRACERR,
                    help="drop fits with tau_err/tau >= this (default %.2f)"
                         % DEFAULT_MAX_FRACERR)
    ap.add_argument("--no-plots", action="store_true",
                    help="print the averaged-rate table but skip saving figures")
    args = ap.parse_args()

    titles = [t for t in CONFIGS if t.startswith(args.prefix)]
    if not titles:
        ap.error("no configs with title prefix %r; use run_exposures.py --list"
                 % args.prefix)

    print("Selecting %d config(s) with prefix %r:" % (len(titles), args.prefix))
    points = load_points(titles, args.indir,
                         args.min_tau, args.min_fracerr, args.max_fracerr)
    if not points:
        ap.error("no usable data found under %s" % args.indir)

    averaged = {temp: average_by_setpoint(rows) for temp, rows in points.items()}
    print_table(averaged)

    if not args.no_plots:
        make_plot(averaged, args.plotdir, args.name)


if __name__ == "__main__":
    main()
