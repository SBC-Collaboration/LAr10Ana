"""Write a dev bubble.sbc for one run by reusing EventDealer's production ProcessSingleRun.

ProcessSingleRun already loops a run's events, loads the images, runs BubbleFinder, and writes
bubble.sbc to the specified output directory. validate.py skips any run with no bubble.sbc.

Make changes to BubbleFinder and re-run this script along with validate.py to
see the effect on one run

    python write_recon.py [run] [out_root] [unpacked_root]

Pass 'all' as the run to regenerate every unique run in the handscan CSV. This is
the long job, usually run on the grid after a single run has been sanity-checked:

    python write_recon.py all [out_root] [unpacked_root] [csv_path]
"""
import csv
import os
import sys

# Weird workaround
_GRID = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "grid_jobs"))
if _GRID not in sys.path:
    sys.path.insert(0, _GRID)

from EventDealer import ProcessSingleRun

HERE = os.path.dirname(os.path.abspath(__file__))  # write recon here, never into shared /exp data

def write_recon(run, out_root, unpacked_root):
    """Run BubbleFinder over `run` -> <out_root>/<run>/bubble.sbc, then return that path."""
    # this is a dev wrapper that deletes/rewrites bubble.sbc. never let it touch shared /exp recon
    if os.path.realpath(out_root).startswith("/exp/"):  # realpath resolves symlinks into /exp
        raise SystemExit(f"refusing to write recon under shared /exp: {out_root}")
    out_dir = os.path.join(out_root, run)
    os.makedirs(out_dir, exist_ok=True)  # ProcessSingleRun only mkdirs the leaf, not parents
    bubble_sbc = os.path.join(out_dir, "bubble.sbc")
    if os.path.exists(bubble_sbc):
        os.remove(bubble_sbc)  # the sbc Writer appends; a stale file would double the rows

    ProcessSingleRun(
        rundir=os.path.join(unpacked_root, run),
        recondir=out_dir,
        process_list=["bubble"],
    )

    print(f"\nwrote {bubble_sbc}")
    return bubble_sbc


def write_recon_all(out_root, unpacked_root, csv_path):
    """Regenerate bubble.sbc for every unique run in the handscan CSV. Long job.

    One run failing (missing images, bad data) logs and is skipped so the batch
    finishes; failures are summarized at the end.
    """
    with open(csv_path) as fh:
        runs = sorted({row["run"] for row in csv.DictReader(fh)})

    print(f"{len(runs)} runs from {csv_path}\n")
    failed = []
    for i, run in enumerate(runs, 1):
        print(f"[{i}/{len(runs)}] {run}")
        try:
            write_recon(run, out_root, unpacked_root)
        except Exception as exc:  # one bad run shouldn't kill the whole batch
            print(f"[fail] {run}: {exc}")
            failed.append(run)
            # drop any partial bubble.sbc so validate.py won't trust a half-written file
            partial = os.path.join(out_root, run, "bubble.sbc")
            if os.path.exists(partial):
                os.remove(partial)

    print(f"\ndone: {len(runs) - len(failed)}/{len(runs)} ok")
    if failed:
        print(f"failed: {', '.join(failed)}")


if __name__ == "__main__":
    # writes recon next to the CSVs (never into shared /exp); raw images read from unpacked
    run = sys.argv[1] if len(sys.argv) > 1 else "20260213_4"
    out_root = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "dev-output")
    unpacked_root = sys.argv[3] if len(sys.argv) > 3 else "/exp/e961/data/SBC-25-unpacked"
    if run == "all":
        csv_path = sys.argv[4] if len(sys.argv) > 4 else os.path.join(HERE, "all_handscans.csv")
        write_recon_all(out_root, unpacked_root, csv_path)
    else:
        write_recon(run, out_root, unpacked_root)
