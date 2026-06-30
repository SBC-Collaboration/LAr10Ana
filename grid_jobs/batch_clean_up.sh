#!/bin/bash
# Batch clean up script
# find all output folders (with or without tag) and runs clean_up.sh with appropriate tags

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLEAN_UP_SCRIPT="${SCRIPT_DIR}/clean_up.sh"

OUTPUT_DIR="/pnfs/coupp/scratch/users/${USER}"
if [ "$USER" = "coupppro" ]; then
    OUTPUT_DIR="/pnfs/coupp/scratch/coupppro"
fi

echo -e "\nStarting batch clean up of grid job outputs in $OUTPUT_DIR"

# Find all grid_output folders
for FOLDER in "$OUTPUT_DIR"/grid_output*; do
    [ -d "$FOLDER" ] || continue
    folder_name=$(basename "$FOLDER")

    # Extract tag from folder name (grid_output or grid_output_<tag>)
    tag="${folder_name#grid_output_}"
    if [ "$tag" = "$folder_name" ]; then
        tag=""
    fi

    echo -e "\nCleaning outputs files for tag '$tag'"

    # Delete files directly under the grid_output folder (corrupted folders)
    find "$FOLDER" -maxdepth 1 -type f -print -delete

    # Clean up
    bash "$CLEAN_UP_SCRIPT" -t "$tag"
done
