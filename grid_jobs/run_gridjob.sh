#!/bin/bash

# Quit if error happens, but don't kill the entire terminal in jupyter.
(
set -e

# Parse options
PRODUCTION_MODE=false
VERBOSE=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --pro|--production|-p) PRODUCTION_MODE=true; shift ;;
        --verbose|-v) VERBOSE=true; shift ;;
        *) RUN_ID=$1; shift ;;
    esac
done

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
OUT_DIR="${DEST_DIR}/grid_output"
# DIR to save a list of jobs. Cannot be on PNFS because it doesn't support append
LIST_DIR="${HOME}/.cache/sbc_job_list.csv"
mkdir -p "$TEMP_DIR"
mkdir -p "$OUT_DIR"
if [ $VERBOSE = true ]; then
    echo "Preparing files for grid job submission for run ${RUN_ID}."
fi

# Copy run tar file over
rsync -azh --chmod=777 "${DATA_DIR}/${RUN_ID}.tar" "${TEMP_DIR}/"
# Tar LAr10ana into a tarball
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LAR10ANA_DIR="${SCRIPT_DIR}/.."
cd "${LAR10ANA_DIR}"
VERSION_FILE="version.txt"
git describe --tags --always >${VERSION_FILE}
TARBALL="LAr10ana.tar"
tar --mtime='1970-01-01 00:00:00' --sort=name -cf $TARBALL --exclude='*.pyc' *.py *.sh ana grid_jobs ${VERSION_FILE}
rm ${VERSION_FILE}
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
# RAM: 3x size + 2GB, minimum 2GB, max 128GB
RAM_GB=$((TAR_SIZE_GB * 3 + 2))
RAM_GB=$((RAM_GB < 2 ? 2 : RAM_GB))
RAM_GB=$((RAM_GB > 128 ? 128 : RAM_GB))
# Run time: 1h/5GB x size + 2h, minimum 2h, maximum 24h
RUN_TIME=$((TAR_SIZE_GB / 5 + 2))
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
	"${TEMP_DIR}/${RUN_ID}.tar" \
    "${OUT_DIR}" \
)
if [ $VERBOSE = true ]; then
    echo -e "\nJob submission output:"
    echo "$output"
fi

# Retrieve Job ID from output, and save to list
JOB_ID=$(echo "$output" | grep -oP '\d+\.\d+@\S+\.fnal\.gov')
echo "$(date '+%Y-%m-%d %H:%M:%S'), ${RUN_ID}, ${JOB_ID}" >> "${LIST_DIR}"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Job ${JOB_ID} for run ${RUN_ID} ($TAR_SIZE_GB GB) successfully submitted."
rm $TARBALL
rm -f "LAr10ana.tar*.tbz"  # remove the tbz files
)

if [[ $? -ne 0 ]]; then
  echo "An error occurred. Quitting."
  rm -f "LAr10ana.tar"
  rm -f "LAr10ana.tar*.tbz"
fi
