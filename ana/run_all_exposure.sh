#!/bin/bash
#
# Run the SBC-25 exposure analysis (run_exposures.py) over every run
# configuration in configs.py, in parallel.
#
# Outputs (defaults baked into run_exposures.py):
#   text  -> /exp/e961/data/users/gputnam/exposures-refactor
#   plots -> /exp/e961/app/users/gputnam/LAr10Ana/plots_06_30_26
#
# Usage:
#   ./run_all.sh [PARALLEL]
# where PARALLEL is the number of concurrent configs (default: nproc-1, leaving
# one core of headroom on this shared interactive node).
#
set -euo pipefail

source /cvmfs/coupp.opensciencegrid.org/LAr10Ana/miniforge3/etc/profile.d/conda.sh
conda activate env

cd "$(dirname "$0")"

NPROC="$(nproc)"
PARALLEL="${1:-$(( NPROC > 1 ? NPROC - 1 : 1 ))}"
LOGDIR="${LOGDIR:-/exp/e961/app/users/gputnam/LAr10Ana/ana/_logs}"
mkdir -p "$LOGDIR"

echo "nproc=$NPROC  parallel=$PARALLEL  logs=$LOGDIR"

# One config per line (titles never contain newlines); xargs -d'\n' keeps each
# whole title (which may contain spaces, '#', quotes, commas) as a single arg.
# -I {} execs directly (no shell), so special characters in titles are safe.
python run_exposures.py --list \
  | xargs -d '\n' -P "$PARALLEL" -I {} \
      bash -c 'python run_exposures.py --config "$1" > "$2/$(echo "$1" | tr " /#,\x27" "_____").log" 2>&1 && echo "OK   $1" || echo "FAIL $1"' _ {} "$LOGDIR"

echo "All configurations processed. Per-config logs in $LOGDIR"
