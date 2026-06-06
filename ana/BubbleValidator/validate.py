"""Pull bubble.sbc, run each nbub algo over a parameter grid, compare against hand scans."""

import functools
import itertools
import os
import sys

import numpy as np
import pandas as pd
from sbcbinaryformat import Streamer

HERE = os.path.dirname(os.path.abspath(__file__))

from nbub import (
    clusters_dbscan,
    clusters_greedy,
    clusters_hdbscan,
    clusters_nearest,
    clusters_radius,
    clusters_significance,
    clusters_track,
    nbub,
)

# every parameter combination in a grid is validated against hand scans
GREEDY_GRID = {
    "max_dist": [10, 15, 20, 25, 30, 35, 40, 45, 50],
    "min_frames": [1, 2, 3, 4, 5],
}
NEAREST_GRID = {
    "max_dist": [10, 15, 20, 25, 30, 35, 40, 45, 50],
    "min_frames": [1, 2, 3, 4, 5],
}
TRACK_GRID = {
    "max_dist": [10, 15, 20, 25, 30, 35, 40, 45, 50],
    "min_frames": [1, 2, 3, 4, 5],
}
TRACK_SIGNIFICANCE_GRID = {
    "max_dist": [25, 30, 35, 40, 45],
    "min_frames": [1, 2, 3],
    "min_sig": [0.0, 0.2, 0.3, 0.4, 0.5],
}
# sklearn density clustering on (x, y) detections
DBSCAN_GRID = {
    "max_dist": [15, 20, 25, 30, 35],
    "min_frames": [1, 2, 3],
}
HDBSCAN_GRID = {
    "max_dist": [15, 20, 25, 30, 35],
    "min_frames": [1, 2, 3],
}
# radius-aware merge; here max_dist is slack added on top of the two bubble radii
RADIUS_GRID = {
    "max_dist": [0, 5, 10, 15, 20, 25, 30],
    "min_frames": [1, 2, 3, 4, 5],
}


def events_from_recon(sbc_path):
    """bubble.sbc -> {ev: per-event dict in nbub's expected shape}"""
    cols = Streamer(sbc_path).to_dict()
    ev_col = cols["ev"]
    out = {}
    for ev in np.unique(ev_col):
        mask = ev_col == ev
        out[int(ev)] = {k: col[mask] for k, col in cols.items()}
    return out


def run_grid(algo, grid, events, scans):
    """Run one algo over its grid; return result rows tagged with algo name + params."""
    keys = list(grid)
    rows = []
    # for each combo of params run the algo
    for combo in itertools.product(*grid.values()):
        params = dict(zip(keys, combo))
        label = ",".join(f"{k}={v}" for k, v in params.items())
        for _, handscan in scans.iterrows():
            if handscan["run"] not in events:
                continue
            event = events[handscan["run"]].get(int(handscan["ev"]))
            # pass partially applied algo with params to nbub
            count = (
                nbub(event, functools.partial(algo, **params))
                if event is not None
                else 0
            )
            rows.append(
                {
                    **handscan.to_dict(),
                    "algo": algo.__name__,
                    "params": label,
                    "reco_nbub": count,
                    "agree": count == int(handscan["scan_nbub"]),
                }
            )
    return rows


def validate(recon_root, csv_path, out_csv):
    handscans = pd.read_csv(csv_path, dtype={"run": str})

    # load events
    events = {}
    for run in handscans["run"].unique():
        path = os.path.join(recon_root, run, "bubble.sbc")
        if os.path.isfile(path):
            events[run] = events_from_recon(path)
        else:
            print(f"[skip] missing {path}")

    # run each algo + param combo, collect rows of results
    rows = []
    rows += run_grid(clusters_greedy, GREEDY_GRID, events, handscans)
    rows += run_grid(clusters_nearest, NEAREST_GRID, events, handscans)
    rows += run_grid(clusters_track, TRACK_GRID, events, handscans)
    rows += run_grid(clusters_significance, TRACK_SIGNIFICANCE_GRID, events, handscans)
    rows += run_grid(clusters_dbscan, DBSCAN_GRID, events, handscans)
    rows += run_grid(clusters_hdbscan, HDBSCAN_GRID, events, handscans)
    rows += run_grid(clusters_radius, RADIUS_GRID, events, handscans)

    res = pd.DataFrame(rows)
    res.to_csv(out_csv, index=False)

    # one row per (algo, param set): how many events scored, and % that matched the handscan
    overall = res.groupby(["algo", "params"])["agree"].agg(n="size", agree_pct="mean")
    overall["agree_pct"] = (overall["agree_pct"] * 100).round(1)
    # same accuracy, but broken out into a column per handscan bubble count (1, 2, 3, ...)
    by_nbub = (
        res.groupby(["algo", "params", "scan_nbub"])["agree"]
        .mean()
        .mul(100)
        .round(1)
        .unstack("scan_nbub")
    )

    print(f"{len(res)} rows -> {out_csv}\n")
    print("accuracy per algo + param set:")
    print(overall.to_string())
    print("\naccuracy (%) per algo + param set, by handscan nbub:")
    print(by_nbub.to_string())
    return res


if __name__ == "__main__":
    # default reads our local recon; pass /exp/e961/data/SBC-25-recon/dev-output for the real one
    recon_root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "dev-output")
    validate(recon_root,
             os.path.join(HERE, "all_handscans.csv"),
             os.path.join(HERE, "validation_results.csv"))
