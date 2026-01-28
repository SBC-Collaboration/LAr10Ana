#!/bin/bash
(
set -e

# Directory where data are stored
DATA_DIR="/exp/e961/data/SBC-25-daqdata"
# Directory where finished job outputs are copied to
RECON_DIR="/exp/e961/data/SBC-25-recon/dev-output"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# File containing all submitted jobs
JOBS_LIST="${HOME}/.cache/sbc_job_list.csv"
echo "Starting batch submission of grid jobs..."

# Parse command line options
FORCE_RERUN=false
VERBOSE=false
PRODUCTION_MODE=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --force|--force-rerun|-f) FORCE_RERUN=true; shift ;;
        --pro|--production|-p) PRODUCTION_MODE=true; shift ;;
        --verbose|-v) VERBOSE=true; shift ;;
        *) RUN_ID=$1; shift ;;
    esac
done

if [ "$FORCE_RERUN" = true ]; then
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
for ((i=${total}-1; i>=0; i--)); do
    tar_file="${tar_files[$i]}"
    # Extract run_id from filename and submit job
    run_id=$(basename "$tar_file" .tar)
    count=$((count + 1))

    EXISTING_VERSION=""
    EXISTING_TAG=""
    CASE=0
    
    # Check if output directory exists and compare versions
    run_output_dir="${RECON_DIR}/${run_id}"
    if [ "$FORCE_RERUN" = false ]; then 
        # Check if job is already submitted
        if [ -f "$JOBS_LIST" ] && grep -q ", ${run_id}, " "$JOBS_LIST"; then
            CASE=1 
        # if no current job, check version of existing output
        elif [ -d "$run_output_dir" ]; then
            version_file="${run_output_dir}/version.txt"
            if [ -f "$version_file" ]; then
                EXISTING_VERSION=$(cat "$version_file")
                EXISTING_TAG="${EXISTING_VERSION%%-*}"
                if [[ "$EXISTING_TAG" == "$CURRENT_TAG" ]]; then
                    CASE=2
                elif [[ "$(printf '%s\n' "$CURRENT_TAG" "$EXISTING_TAG" | sort -V | head -n1)" == "$CURRENT_TAG" ]]; then
                    CASE=3
                else
                    CASE=4  # current version is newer
                fi
            else
                CASE=5  # no version.txt
            fi
        else
            CASE=6  # no existing output
        fi
    else
        CASE=7  # force rerun
    fi

    # Log if verbose or not skipping
    SKIPPING=false
    if [ "$CASE" -eq 1 ] || [ "$CASE" -eq 2 ] || [ "$CASE" -eq 3 ]; then
        SKIPPING=true
    fi
    if [ "$VERBOSE" = true ] || [ "$SKIPPING" = false ]; then
        echo -e "\n========== $(date '+%Y-%m-%d %H:%M:%S') =========="
        echo "Processing ${count}/${total}: Run ${run_id}"
        if [ -n "$EXISTING_VERSION" ] && [ "$VERBOSE" = true ]; then
            echo "Found existing output with version: ${EXISTING_VERSION} (tag: ${EXISTING_TAG})"
        fi
    fi

    case $CASE in
        1) 
            if [ "$VERBOSE" = true ]; then
                echo "Skipping ${run_id}: job already submitted (found in job list)"
            fi
            ;;
        2) 
            if [ "$VERBOSE" = true ]; then
                echo "Skipping ${run_id}: existing version is same"
            fi
            ;;
        3) 
            if [ "$VERBOSE" = true ]; then
                echo "Skipping ${run_id}: existing version (${EXISTING_TAG}) is newer than current (${CURRENT_TAG})"
            fi
            ;;
        4) 
            echo "Proceeding: current version (${CURRENT_VERSION}) is newer than existing (${EXISTING_VERSION})"
            ;;
        5) 
            echo "Proceeding: no version.txt found in existing output"
            ;;
        6) 
            echo "Proceeding: no existing output found"
            ;;
        7) 
            echo "Proceeding: Force rerun enabled"
            ;;
        *) echo "Unknown case $CASE"
            ;;
    esac
    
    if [ "$SKIPPING" = false ]; then
        if [ "$VERBOSE" = true ]; then
            VERBOSE_ARG="--verbose"
        else
            VERBOSE_ARG=""
        fi

        if [ "$PRODUCTION_MODE" = true ]; then
            ROLE_ARG="--production"
        else
            ROLE_ARG=""
        fi

        bash "${SCRIPT_DIR}/run_gridjob.sh" ${VERBOSE_ARG} ${ROLE_ARG} "$run_id"
    fi
done

echo "Batch submission of ${total} jobs complete."
)

if [[ $? -ne 0 ]]; then
    echo "An error occurred during batch submission. Quitting."
fi