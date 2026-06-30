#!/bin/bash

# Quit if error happens, but don't kill the entire terminal in jupyter.
(
set -e

# Parse options
PRODUCTION_MODE=false
VERBOSE=false
TAG=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --pro|--production|-p) PRODUCTION_MODE=true; shift ;;
        --verbose|-v) VERBOSE=true; shift ;;
        --tag|-t) 
            TAG="${2:-}"
            [ -n "$TAG" ] || { echo "Error: --tag requires a value (e.g. --tag v0.1.0)"; exit 2; }
            shift 2 ;;
        *) RUN_ID=$1; shift ;;
    esac
done

if [ -z "${RUN_ID:-}" ]; then
    echo "Usage: $0 [--tag vX.Y.Z] [--production] [--verbose] RUN_ID"
    exit 2
fi

ROLE_ARG=""
if [ "$PRODUCTION_MODE" = true ]; then
    ROLE_ARG="--role production"
fi

# setup jobsub environment
JOBSUB_LITE_SH="/etc/profile.d/jobsub_lite.sh"
if [ -f "$JOBSUB_LITE_SH" ]; then
    source "$JOBSUB_LITE_SH"
fi
export HTGETTOKENOPTS="--credkey=coupppro/managedtokens/fifeutilgpvm01.fnal.gov"

DEST_DIR="/pnfs/coupp/scratch/users/${USER}"
if [ "$USER" = "coupppro" ]; then
    DEST_DIR="/pnfs/coupp/scratch/coupppro"
fi
DATA_DIR="/exp/e961/data/SBC-25-daqdata"
# Directory where run data will be copied to for this grid job
TEMP_DIR="${DEST_DIR}/temp_data"
# Directory where job output will be saved to
OUT_SUFFIX=""
SAFE_TAG=""
if [ -n "$TAG" ]; then
    SAFE_TAG="$(echo "$TAG" | tr '/ ' '_')"
    OUT_SUFFIX="_${SAFE_TAG}"
fi
OUT_DIR="${DEST_DIR}/grid_output${OUT_SUFFIX}"
# DIR to save a list of jobs. Cannot be on PNFS because it doesn't support append
LIST_FILE="${HOME}/.cache/sbc/jobs_list.csv"
if [ -n "$TAG" ]; then
    LIST_FILE="${HOME}/.cache/sbc/jobs_list_${SAFE_TAG}.csv"
fi
mkdir -p "$(dirname "$LIST_FILE")"
mkdir -p "$TEMP_DIR"
mkdir -p "$OUT_DIR"
if [ $VERBOSE = true ]; then
    echo "Preparing files for grid job submission for run ${RUN_ID}."
fi

SRC="${DATA_DIR}/${RUN_ID}.tar"
DST="${TEMP_DIR}/${RUN_ID}.tar"

# Sometimes the PNFS can be corrupted, so do a quick read to check if IO is working.
if [[ -f "$DST" ]] && ! head -c 1M "$DST" > /dev/null 2>&1; then
  echo "BAD PNFS copy, removing: $DST"
  rm -f "$DST"
fi

# Copy run tar file over
rsync -azh --chmod=777 "${SRC}" "${DST}"
# Tar LAr10ana into a tarball (at the requested tag)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LAR10ANA_DIR="${SCRIPT_DIR}/.."
VERSION_FILE="version.txt"
TARBALL="LAr10ana.tar"

cleanup() {
    # Always clean tar/tbz
    rm -f "${LAR10ANA_DIR}/${TARBALL}" "${LAR10ANA_DIR}/LAr10ana.tar" "${LAR10ANA_DIR}/LAr10ana.tar"*.tbz* >/dev/null 2>&1 || true
    # Clean worktree if it exists
    if [ -n "${WT_DIR:-}" ]; then
        git -C "${LAR10ANA_DIR}" worktree remove --force "${WT_DIR}" >/dev/null 2>&1 || true
        rm -rf "${WT_DIR}" >/dev/null 2>&1 || true
    fi
}
trap cleanup EXIT

if [ -n "$TAG" ]; then
    # Update to latest tags
    git -C "${LAR10ANA_DIR}" fetch --tags --quiet || true
    # Create a temporary worktree folder
    WT_DIR="$(mktemp -d -t lar10ana_wt_XXXXXX)"
    # Add the worktree at the requested tag
    git -C "${LAR10ANA_DIR}" worktree add --detach "${WT_DIR}" "refs/tags/${TAG}" >/dev/null
    ( cd "${WT_DIR}"
      # Inside the worktree, generate version file
      git describe --tags --always > "${VERSION_FILE}"
      # Tar the required files
      tar --mtime='1970-01-01 00:00:00' --sort=name -cf "${LAR10ANA_DIR}/${TARBALL}" --exclude='*.pyc' *.py *.sh ana grid_jobs "${VERSION_FILE}"
      rm -f "${VERSION_FILE}"
    )
else
    cd "${LAR10ANA_DIR}"
    # Generate version file and tar required files
    git describe --tags --always >${VERSION_FILE}
    tar --mtime='1970-01-01 00:00:00' --sort=name -cf $TARBALL --exclude='*.pyc' *.py *.sh ana grid_jobs ${VERSION_FILE}
    rm ${VERSION_FILE}
fi
if [ $VERBOSE = true ]; then
    echo "Data copied over. LAr10ana is tarred. Ready for job submission."
fi

# Calculate disk and memory usage by tar file size
TAR_SIZE_BYTES=$(stat -f%z "${DATA_DIR}/${RUN_ID}.tar" 2>/dev/null || stat -c%s "${DATA_DIR}/${RUN_ID}.tar")
TAR_SIZE_GB=$((TAR_SIZE_BYTES / 1024 / 1024 / 1024))
# Disk: 3x size + 5GB, min 5GB, max 500GB
DISK_GB=$((TAR_SIZE_GB * 3 + 5))
DISK_GB=$((DISK_GB < 5 ? 5 : DISK_GB))
DISK_GB=$((DISK_GB > 500 ? 500 : DISK_GB))
# RAM: 2x size + 2GB, minimum 2GB, max 16GB
RAM_GB=$((TAR_SIZE_GB * 2 + 2))
RAM_GB=$((RAM_GB < 2 ? 2 : RAM_GB))
RAM_GB=$((RAM_GB > 16 ? 16 : RAM_GB))
# Run time: 1h/2GB x (size + 1) + 2h, minimum 2h, maximum 24h
RUN_TIME=$(( (TAR_SIZE_GB + 1) / 2 + 2))
RUN_TIME=$((RUN_TIME < 2 ? 2 : RUN_TIME))
RUN_TIME=$((RUN_TIME > 24 ? 24 : RUN_TIME))
if [ $VERBOSE = true ]; then
    echo "Tar size: ${TAR_SIZE_GB}GB, requesting Disk: ${DISK_GB}GB, RAM: ${RAM_GB}GB, Run Time ${RUN_TIME}h"
fi

# Submit job and pipe output
output=$( \
  jobsub_submit --disk=${DISK_GB}GB --expected-lifetime=${RUN_TIME}h --memory=${RAM_GB}GB -G coupp \
    --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE,DEDICATED \
	--tar_file_name dropbox:///${LAR10ANA_DIR}/${TARBALL} \
	-N 1 ${ROLE_ARG} \
	file://${SCRIPT_DIR}/gridjob.sh \
	"${DST}" \
    "${OUT_DIR}" \
)
if [ $VERBOSE = true ]; then
    echo -e "\nJob submission output:"
    echo "$output"
fi

# Retrieve Job ID from output, and save to list
JOB_ID=$(echo "$output" | grep -oP '\d+\.\d+@\S+\.fnal\.gov')
echo "$(date '+%Y-%m-%d %H:%M:%S'), ${RUN_ID}, ${JOB_ID}" >> "${LIST_FILE}"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Job ${JOB_ID} for run ${RUN_ID} ($new_version, $TAR_SIZE_GB GB) successfully submitted."
)
