#!/usr/bin/env python
"""Standalone SBC-25 (LAr10) exposure analysis.

Refactor of ExposureAnalysis.ipynb + mix_exposures.py into a single script that,
for one run configuration:

  * reads the pre-processed recon outputs (exposure.sbc + event.sbc) instead of
    re-running the manual pressure/livetime computation on raw DAQ data,
  * fits an exponential to the per-expansion livetime distribution at each
    pressure setpoint to extract a lifetime tau,
  * directly produces the final "mix" result -- the nominal (coarse-binning) fit
    at each setpoint, falling back to the zoom (fine-binning) fit when the
    nominal lifetime is smaller than the nominal bin width -- with no
    intermediate nominal/zoom files,
  * makes one livetime-fit plot per setpoint,
  * appends physics columns from the SeitzModel library: number of events in the
    fit window, Seitz threshold Q, Eion, rion, and rho_l (pho_l).

The fit machinery (dofitexp) and the pressure-setpoint snapping reproduce the
notebook exactly, so for runs present in recon the output reproduces the legacy
results in /exp/e961/data/users/gputnam/exposures_1-29 (see --help / README).

Fit parameters match the exposures_1-29 reference:
    nominal : THI=180, TLO=0,  bins = linspace(0, 300, 16)  (20 s wide)
    zoom    : THI=40,  TLO=6,  bins = linspace(0,  40, 21)  ( 2 s wide)
    mix rule: use the nominal (large-bin) fit; if its lifetime is smaller than
              the nominal bin width (20 s) AND the zoom window holds >=1 event,
              use the zoom (small-bin) fit instead.

Usage:
    python run_exposures.py --config "Cold Cs 1/9"
    python run_exposures.py --all
"""

import argparse
import os
import sys
import warnings

import numpy as np

import matplotlib
matplotlib.use("Agg")          # headless: no display needed for batch plotting
import matplotlib.pyplot as plt

from scipy.optimize import minimize, minimize_scalar

from configs import CONFIGS, sanitize

warnings.filterwarnings("ignore", category=RuntimeWarning)

# --- SeitzModel (use the CVMFS copy: it has the .E and .Rl attributes) ---------
sys.path.append("/cvmfs/coupp.opensciencegrid.org/LAr10Ana/SeitzModel/REFPROP")
import SeitzModel as sm

from sbcbinaryformat import Streamer

# --- fixed analysis constants (match the exposures_1-29 reference) -------------
REQUIRE_QUIET = False
LIVETIME_CUT = -1
FLUID = "argon"
BAR_TO_PSIA = 14.5038

# Fit windows / binning per flavor.
NOMINAL = dict(thi=180, tlo=0, bins=np.linspace(0, 300, 16))
ZOOM = dict(thi=40, tlo=6, bins=np.linspace(0, 40, 21))

# Mixing rule: at each setpoint do the nominal (large-bin) fit; if its fitted
# lifetime is smaller than the nominal bin width -- the decay happens within a
# single large bin and cannot be resolved -- use the zoom (small-bin) fit for
# that setpoint instead.
NOMINAL_BINWIDTH = NOMINAL["bins"][1] - NOMINAL["bins"][0]   # 20 s

# Operating temperature rule (kelvin), keyed on run date.
TEMP_DATE_CUT = 20260130        # 1/30/2026
TEMP_BEFORE_K = 116.7
TEMP_AFTER_K = 119.6

DEFAULT_RECONDIR = "/exp/e961/data/SBC-25-recon/v0.2.0"
DEFAULT_OUTDIR = "/exp/e961/data/users/gputnam/exposures-refactor"
DEFAULT_PLOTDIR = "/exp/e961/app/users/gputnam/LAr10Ana/plots_06_30_26"

OUT_HEADER = (
    "Pressure [bara]\tLifetime [s]\tLifetime Error [s]\t"
    "Exponential Fit 2xNLL\tN.d.o.f.\tTime Cut High [s]\tTime Cut Low [s]\t"
    "N events\tSeitz Threshold [keV]\tEion [keV]\trion [nm]\tpho_l [g/cc]"
)


# ----------------------------------------------------------------------------- #
# Helpers reproduced from ExposureAnalysis.ipynb
# ----------------------------------------------------------------------------- #
def snap(v, step):
    """Round to the nearest multiple of step (pressure-setpoint binning)."""
    return (v / step).round() * step


def dofitexp(X, N, tau0=50, tlo=0, thi=180):
    """Exponential fit to a binned livetime histogram (verbatim from the notebook).

    Returns (tau, tau_err, 2*nll_min, ndof, evalf). Poisson likelihood-ratio
    cost, profiled 1-sigma error. Raises RuntimeError on a failed minimize.
    """
    tau_width = X[1] - X[0]

    N = N[(X < thi) & (X > tlo)]
    X = X[(X < thi) & (X > tlo)]

    def fitexp(X, tau):
        tlo_norm = max(LIVETIME_CUT, tlo)
        norm = tau * (np.exp(-tlo_norm / tau) - np.exp(-thi / tau))
        return (N.sum() * tau_width / norm) * np.exp(-X / tau)

    def nll(p):
        mu = fitexp(X, *p)
        if np.any(mu < 0):
            return np.inf
        return 2 * (np.sum(mu - N) + np.sum(N[N > 0] * np.log(N[N > 0] / mu[N > 0])))

    # Robust 1-D bounded minimization. The previous L-BFGS-B-from-fixed-tau0
    # approach stalled near the start point for ~1/3 of setpoints (finite-
    # difference gradient + line search terminating early), returning an
    # inflated tau and a visibly poor fit. minimize_scalar's bounded (Brent)
    # search needs no seed and finds the true minimum for every setpoint.
    result = minimize_scalar(lambda t: nll([t]), bounds=(0.01, 1e4),
                             method="bounded")
    tau = result.x
    nll_min = nll([tau])

    dtau = (np.linspace(0, 2, 1000) * tau)[1:]
    profile = np.array([nll([dt]) for dt in dtau]) - nll_min
    tau_1sigma = dtau[np.argmin(np.abs(profile - 1))]
    tau_err = np.abs(tau - tau_1sigma)

    return tau, tau_err, nll_min, X.size, lambda X: fitexp(X, tau)


# ----------------------------------------------------------------------------- #
# Data loading from recon outputs
# ----------------------------------------------------------------------------- #
def load_config(runs, recondir):
    """Load per-event arrays for a configuration from recon exposure/event files.

    For each run we inner-join exposure.sbc (PT2121 pressure/livetime, quiet flag)
    with event.sbc (ev_livetime, ev_exit_code) on the per-event index ``ev``.

    Returns (pressures, livetimes, times_ms, code, quiet, present, missing).
    """
    P, L, T, C, Q = [], [], [], [], []
    present, missing = [], []
    for r in runs:
        epath = os.path.join(recondir, r, "exposure.sbc")
        vpath = os.path.join(recondir, r, "event.sbc")
        if not (os.path.exists(epath) and os.path.exists(vpath)):
            missing.append(r)
            continue
        e = Streamer(epath).to_dict()
        v = Streamer(vpath).to_dict()

        # Map event.sbc rows by ev index, then align to exposure.sbc rows.
        ev_v = np.asarray(v["ev"], int)
        vlt = {ev: lt for ev, lt in zip(ev_v, np.asarray(v["ev_livetime"], float))}
        vcd = {ev: cd for ev, cd in zip(ev_v, np.asarray(v["ev_exit_code"], float))}

        ev_e = np.asarray(e["ev"], int)
        pe = np.asarray(e["PT2121_pressure"], float)
        le = np.asarray(e["PT2121_livetime"], float)
        qe = np.asarray(e["is_quiet_mode"], bool)
        for j, ev in enumerate(ev_e):
            if ev not in vlt:
                continue
            P.append(pe[j])
            L.append(le[j])
            Q.append(qe[j])
            T.append(vlt[ev])
            C.append(vcd[ev])
        present.append(r)

    if not P:
        return (np.array([]),) * 5 + (present, missing)
    return (np.array(P), np.array(L), np.array(T), np.array(C),
            np.array(Q, dtype=bool), present, missing)


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


def seitz_values(pressure_bara, temp_K):
    """(Q, Eion, rion, rho_l) from SeitzModel for one pressure / temperature."""
    try:
        s = sm.SeitzModel(pressure_bara * BAR_TO_PSIA, temp_K - 273.15, FLUID)
        return float(s.Q), float(s.E), float(s.Rl), float(s.Rho_l)
    except Exception as exc:               # pragma: no cover - defensive
        print("  WARNING: SeitzModel failed at P=%.3f bara: %r" % (pressure_bara, exc))
        return np.nan, np.nan, np.nan, np.nan


# ----------------------------------------------------------------------------- #
# Per-configuration processing
# ----------------------------------------------------------------------------- #
def _fit_flavor(flavor, lt):
    """Histogram livetimes and fit one flavor (NOMINAL or ZOOM).

    Returns a dict with the histogram, fit window, n_window, and (when the fit
    succeeds) tau/tau_err/nll/ndof/evalf. ``ok`` is False if dofitexp fails.
    """
    thi, tlo, bins = flavor["thi"], flavor["tlo"], flavor["bins"]
    centers = (bins[:-1] + bins[1:]) / 2
    N, _ = np.histogram(lt, bins=bins)
    res = dict(thi=thi, tlo=tlo, bins=bins, centers=centers, N=N,
               n_window=int(N[(centers < thi) & (centers > tlo)].sum()),
               zoomed=flavor is ZOOM, ok=True)
    try:
        tau, tau_err, nll, ndof, evalf = dofitexp(centers, N, thi=thi, tlo=tlo)
        res.update(tau=tau, tau_err=tau_err, nll=nll, ndof=ndof, evalf=evalf)
    except RuntimeError:
        res["ok"] = False
    return res


def _plot_setpoint(plotdir, stem, title, p, res, lt):
    """One livetime-fit figure for a setpoint, using the chosen flavor's result."""
    fig = plt.figure()
    plt.hist(lt, bins=res["bins"])
    plt.title("%s (%.2f bar)" % (title, p))
    plt.xlabel("Livetime [s]")
    plt.ylabel("Events")
    if res["ok"]:
        if not res["zoomed"]:
            plt.axvline([res["thi"]], color="gray", linestyle=":")
        if res["tlo"] > 0:
            plt.axvline([res["tlo"]], color="gray", linestyle=":")
        plt.text(0.6, 0.95,
                 "Exponential Fit\n$\\tau$: %.1f $\\pm$ %.1f [s]\n"
                 "$-2\\Delta\\ln\\mathcal{L}/n: %.1f / %i$"
                 % (res["tau"], res["tau_err"], res["nll"], res["ndof"]),
                 verticalalignment="top", transform=plt.gca().transAxes)
        plt.plot(res["centers"], res["evalf"](res["centers"]), color="red")
        plt.errorbar(res["centers"], res["N"], np.sqrt(res["N"]),
                     linestyle="none", color="black")
    _savefig(fig, plotdir, stem, p, res["zoomed"])
    plt.close(fig)


def process_config(title, runs, recondir, outdir, plotdir, dosave=True):
    stem = sanitize(title)
    print("=== %s  (%d runs) ===" % (title, len(runs)))

    pressures, livetimes, times, code, quiet, present, missing = load_config(runs, recondir)
    if missing:
        print("  MISSING recon (no exposure.sbc/event.sbc): %s" % ", ".join(missing))
    if pressures.size == 0:
        print("  no events loaded -- skipping")
        return None

    valid = (times > 40e3) & ((code == 0) | (code == 2001))
    if REQUIRE_QUIET:
        valid = valid & quiet
    if LIVETIME_CUT > 0:
        valid = valid & (times > LIVETIME_CUT)

    pset = snap(pressures, 0.25)
    ps = np.unique(pset[valid])
    ps = ps[~np.isnan(ps)]

    temp_K = temperature_K(present if present else runs)
    print("  setpoints: %s" % ", ".join("%.2f" % p for p in ps))
    print("  temperature: %.1f K" % temp_K)

    rows = []
    for p in ps:
        lt = livetimes[valid & (pset == p)]

        # Always fit the nominal (large-bin) flavor; if its lifetime is smaller
        # than the nominal bin width, the large bins cannot resolve the decay --
        # use the zoom (small-bin) fit for this setpoint instead, but only when
        # the zoom fit window actually contains at least one event (otherwise an
        # empty window leaves the zoom fit unconstrained).
        res = _fit_flavor(NOMINAL, lt)
        if res["ok"] and res["tau"] < NOMINAL_BINWIDTH:
            zoom = _fit_flavor(ZOOM, lt)
            if zoom["n_window"] >= 1:
                res = zoom

        if plotdir:
            _plot_setpoint(plotdir, stem, title, p, res, lt)

        if not res["ok"]:
            rows.append([p, -1, -1, -1, -1, -1, -1, res["n_window"],
                         *seitz_values(p, temp_K)])
            continue
        rows.append([p, res["tau"], res["tau_err"], res["nll"], res["ndof"],
                     res["thi"], res["tlo"], res["n_window"], *seitz_values(p, temp_K)])

    out = np.array(rows)
    if dosave:
        os.makedirs(outdir, exist_ok=True)
        fout = os.path.join(outdir, "%s_exposures_mix.txt" % stem)
        np.savetxt(fout, out, header=OUT_HEADER, comments="# ")
        print("  wrote %s" % fout)
    return out


def _savefig(fig, plotdir, stem, p, zoomed):
    os.makedirs(plotdir, exist_ok=True)
    suffix = "_exposure_zoom_P%.2f" if zoomed else "_exposure_P%.2f"
    base = (stem + (suffix % p)).replace(" ", "-").replace(".", "_").replace("/", "_")
    fig.savefig(os.path.join(plotdir, base + ".pdf"))
    fig.savefig(os.path.join(plotdir, base + ".png"))


# ----------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--config", help="run configuration TITLE (see configs.py)")
    g.add_argument("--all", action="store_true", help="process every configuration")
    g.add_argument("--list", action="store_true", help="list configuration titles and exit")
    ap.add_argument("--recondir", default=DEFAULT_RECONDIR)
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    ap.add_argument("--plotdir", default=DEFAULT_PLOTDIR)
    ap.add_argument("--no-plots", action="store_true",
                    help="still write the text file but skip saving plots")
    ap.add_argument("--zoom-tlo", type=float, default=None,
                    help="override the zoom-fit lower time cut (default %.1f s; "
                         "use 0 for no lower limit)" % ZOOM["tlo"])
    args = ap.parse_args()

    if args.zoom_tlo is not None:
        print("Overriding zoom TLO: %.1f -> %.1f s" % (ZOOM["tlo"], args.zoom_tlo))
        ZOOM["tlo"] = args.zoom_tlo

    if args.list:
        for t in CONFIGS:
            print(t)
        return

    if args.all:
        titles = list(CONFIGS)
    else:
        if args.config not in CONFIGS:
            ap.error("unknown config %r; use --list to see valid titles" % args.config)
        titles = [args.config]

    for t in titles:
        process_config(t, CONFIGS[t], args.recondir, args.outdir,
                       None if args.no_plots else args.plotdir,
                       dosave=True)


if __name__ == "__main__":
    main()
