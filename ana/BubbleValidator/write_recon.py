"""Write a dev bubble.sbc for one run by reusing EventDealer's production ProcessSingleRun.

ProcessSingleRun already loops a run's events, loads the images, runs BubbleFinder, and writes
bubble.sbc to the specified output directory. validate.py skips any run with no bubble.sbc.

Make changes to BubbleFinder and re-run this script along with validate.py to
see the effect on one run

    python write_recon.py [run] [out_root] [unpacked_root]

Pass 'all' as the run to regenerate every unique run in the handscan CSV. This is
the long job; runs are independent so it fans out across CPU cores (--nprocs):

    python write_recon.py all [out_root] [unpacked_root] [--nprocs N]

Run with -h/--help for all options.
"""
import argparse
import csv
import os
import sys
from multiprocessing import Pool

# Weird workaround
_GRID = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "grid_jobs"))
if _GRID not in sys.path:
    sys.path.insert(0, _GRID)

from EventDealer import ProcessSingleRun

HERE = os.path.dirname(os.path.abspath(__file__))  # write recon here, never into shared /exp data


def write_recon(run, out_root, unpacked_root):
    """Run BubbleFinder over `run` -> <out_root>/<run>/bubble.sbc, then return that path."""
    # this is a dev wrapper that deletes/rewrites bubble.sbc. never let it touch shared /exp recon
    if os.path.realpath(out_root).startswith("/exp/"):
        raise SystemExit(f"refusing to write recon under shared /exp: {out_root}")
    out_dir = os.path.join(out_root, run)
    os.makedirs(out_dir, exist_ok=True)
    bubble_sbc = os.path.join(out_dir, "bubble.sbc")
    if os.path.exists(bubble_sbc):
        os.remove(bubble_sbc)  # the sbc Writer appends; a stale file would double the rows

    ProcessSingleRun(
        rundir=os.path.join(unpacked_root, run),
        recondir=out_dir,
        process_list=["bubble"],
        # maxevt=3,  # uncomment to cap events per run for a quick dev iteration
    )

    print(f"\nwrote {bubble_sbc}")
    return bubble_sbc


def _run_one(task):
    """Pool worker: process one run. Returns (run, ok, err_or_None).

    Each run owns its own <out_root>/<run>/bubble.sbc, so workers never share the appending
    sbc Writer -- that's why run-level is the safe unit of parallelism."""
    run, out_root, unpacked_root = task
    try:
        write_recon(run, out_root, unpacked_root)
        return run, True, None
    except Exception as exc:  # one bad run shouldn't kill the whole batch
        # drop any partial bubble.sbc so validate.py won't trust a half-written file
        partial = os.path.join(out_root, run, "bubble.sbc")
        if os.path.exists(partial):
            os.remove(partial)
        return run, False, str(exc)


def write_recon_all(out_root, unpacked_root, csv_path, nprocs=None):
    """Regenerate bubble.sbc for every unique run in the handscan CSV. Long job.

    Runs are independent, so they fan out across a multiprocessing.Pool. Default
    uses all cores, capped at the run count. nprocs=1 forces the 1 process. One run
    failing (missing images, bad data) logs and is skipped so the batch finishes.
    """
    run_events = {}
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            run_events.setdefault(row["run"], set()).add(int(row["ev"]))
    # start the runs with the most handscan events before the small ones,
    # so the heaviest job isn't left running alone at the end on a single core.
    runs = sorted(run_events, key=lambda r: len(run_events[r]), reverse=True)
    tasks = [(r, out_root, unpacked_root) for r in runs]

    if nprocs is None:
        nprocs = min(len(runs), os.cpu_count() or 1)
    nprocs = max(1, min(nprocs, len(runs)))

    print(f"{len(runs)} runs from {csv_path} -- {nprocs} worker(s)\n")
    failed = []

    def _report(i, run, ok, err):
        if ok:
            print(f"[{i}/{len(runs)}] {run} ok")
        else:
            print(f"[{i}/{len(runs)}] [fail] {run}: {err}")
            failed.append(run)

    if nprocs == 1:
        for i, task in enumerate(tasks, 1):
            _report(i, *_run_one(task))
    else:
        # imap_unordered: report each run as it finishes rather than waiting for the slowest.
        with Pool(processes=nprocs) as pool:
            for i, (run, ok, err) in enumerate(pool.imap_unordered(_run_one, tasks), 1):
                _report(i, run, ok, err)

    print(f"\ndone: {len(runs) - len(failed)}/{len(runs)} ok")
    if failed:
        print(f"failed: {', '.join(failed)}")


def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Regenerate a dev bubble.sbc for one run, or 'all' runs in the handscan CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # writes recon next to the CSVs (never into shared /exp); raw images read from unpacked
    p.add_argument("run", nargs="?", default="20260213_4",
                   help="run name to reprocess, or 'all' for every run in the CSV")
    p.add_argument("out_root", nargs="?", default=os.path.join(HERE, "dev-output"),
                   help="recon output root (refused under shared /exp)")
    p.add_argument("unpacked_root", nargs="?", default="/exp/e961/data/SBC-25-unpacked",
                   help="root of the unpacked raw images")
    p.add_argument("--csv", default=os.path.join(HERE, "all_handscans.csv"),
                   help="handscan CSV providing the run list for 'all'")
    p.add_argument("--nprocs", type=int, default=None,
                   help="parallel workers for 'all' (default: all cores, capped at run count)")
    return p.parse_args(argv)


if __name__ == "__main__":
    a = _parse_args()
    if a.run == "all":
        write_recon_all(a.out_root, a.unpacked_root, a.csv, nprocs=a.nprocs)
    else:
        write_recon(a.run, a.out_root, a.unpacked_root)
