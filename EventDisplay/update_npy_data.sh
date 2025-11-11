#!/bin/bash

# Script to add new event data to raw_events.npy
set -e

# Resolve absolute path of the cron script
NPY_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# LAr10Ana project root
PROJECT_DIR="$(dirname "$NPY_SCRIPT_DIR")"

echo "$(date '+%Y-%m-%dT%H:%M:%S%z')  Starting update_npy_data.sh"
echo "Npy Script dir: $NPY_SCRIPT_DIR"
echo "Project dir: $PROJECT_DIR"

echo "Sourcing conda environment"
source "$PROJECT_DIR/setup.sh"

echo "Running convert_raw_to_npy_run_by_run.py"
python "$NPY_SCRIPT_DIR/convert_raw_to_npy_run_by_run.py" /exp/e961/data/SBC-25-daqdata "$NPY_SCRIPT_DIR/npy/SBC-25"

echo "Running merge_raw_run_npy.py"
python "$NPY_SCRIPT_DIR/merge_raw_run_npy.py" "$NPY_SCRIPT_DIR/npy/SBC-25/"

echo "$(date '+%Y-%m-%dT%H:%M:%S%z')  Finished update_npy_data.sh"
