#!/bin/bash
(
set -e

# Directory where data are stored
DATA_DIR="/exp/e961/data/SBC-25-daqdata"
# Directory where finished job outputs are copied to
RECON_DIR="/exp/e961/data/SBC-25-recon/dev-output"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo "Starting batch submission of grid jobs..."

# Parse command line options
force_rerun=false
while getopts "f" opt; do
  case $opt in
    f)
      force_rerun=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if [ "$force_rerun" = true ]; then
    echo "Force rerun enabled: Skipping version checks."
fi

# Get current version
CURRENT_VERSION=$(git describe --tags --always)
CURRENT_TAG="${CURRENT_VERSION%%-*}"
echo "Current version: ${CURRENT_VERSION} (tag: ${CURRENT_TAG})"

# Count total tar files
REGEX="${DATA_DIR}"/*.tar
mapfile -t tar_files < <(ls -1 ${REGEX} 2>/dev/null)
total=${#tar_files[@]}

if [ $total -eq 0 ]; then
    echo "No .tar files found in ${DATA_DIR}"
    exit 1
fi
echo "Found ${total} .tar files to process"
count=0

# Loop through all .tar files in DATA_DIR
# in reverse order, since later runs tend to take longer
# for ((i=${total}-1; i>=0; i--)); do
for ((i=0; i<=${total}-1; i++)); do
    tar_file="${tar_files[$i]}"
    # Extract run_id from filename and submit job
    run_id=$(basename "$tar_file" .tar)
    count=$((count + 1))
    echo -e "\n=========================================="
    echo "Processing ${count}/${total}: Run ${run_id}"
    
    # Check if output directory exists and compare versions
    run_output_dir="${RECON_DIR}/${run_id}"
    if [ "$force_rerun" = false ] && [ -d "$run_output_dir" ]; then
        version_file="${run_output_dir}/version.txt"
        if [ -f "$version_file" ]; then
            existing_version=$(cat "$version_file")
            existing_tag="${existing_version%%-*}"
            echo "Found existing output with version: ${existing_version} (tag: ${existing_tag})"
            if [[ "$existing_tag" == "$CURRENT_TAG" ]]; then
                echo "Skipping ${run_id}: existing version is same"
                continue
            elif [[ "$(printf '%s\n' "$CURRENT_TAG" "$existing_tag" | sort -V | head -n1)" == "$CURRENT_TAG" ]]; then
                echo "Skipping ${run_id}: existing version (${existing_tag}) is newer than current (${CURRENT_TAG})"
                continue
            else
                echo "Proceeding: current version (${CURRENT_VERSION}) is newer than existing (${existing_version})"
            fi
        else
            echo "Proceeding: no version.txt found in existing output"
        fi
    elif [ "$force_rerun" = true ]; then
         echo "Proceeding: Force rerun enabled"
    else
        echo "Proceeding: no existing output found"
    fi
    
    "${SCRIPT_DIR}/run_gridjob.sh" "$run_id"
done

echo "Batch submission of ${total} jobs complete."
)

if [[ $? -ne 0 ]]; then
    echo "An error occurred during batch submission. Quitting."
fi