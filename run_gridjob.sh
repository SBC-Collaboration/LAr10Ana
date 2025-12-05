#!/bin/bash

# Quit if error happens, but don't kill the entire terminal in jupyter.
(
set -e

RUN_ID=$1
DATA_DIR="/exp/e961/data/SBC-25-daqdata"
# Directory where run data will be copied to for this grid job
TEMP_DIR="/pnfs/coupp/scratch/users/${USER}/temp_data"
# Directory where job output will be saved to
OUT_DIR="/pnfs/coupp/scratch/users/${USER}/grid_output"
# DIR to save a list of jobs. Cannot be on PNFS because it doesn't support append
LIST_DIR="${HOME}/.cache/sbc_job_list.csv"
mkdir -p "$TEMP_DIR"
mkdir -p "$OUT_DIR"
echo "Preparing files for grid job submission for run ${RUN_ID}."

# Copy run tar file over
rsync -azh --chmod=755 "${DATA_DIR}/${RUN_ID}.tar" "${TEMP_DIR}/"
# Tar LAr10ana into a tarball
TARBALL="LAr10ana_${RUN_ID}.tar"
tar -cf $TARBALL --exclude='*.pyc' *.py ana setup.sh
echo "Data copied over. LAr10ana is tarred. Ready for job submission."

# Calculate disk and memory usage by tar file size
TAR_SIZE_BYTES=$(stat -f%z "${DATA_DIR}/${RUN_ID}.tar" 2>/dev/null || stat -c%s "${DATA_DIR}/${RUN_ID}.tar")
TAR_SIZE_GB=$((TAR_SIZE_BYTES / 1024 / 1024 / 1024))
# Disk: 3x data + 5GB, min 10GB, max 300GB
DISK_GB=$((TAR_SIZE_GB * 2 + 5))
DISK_GB=$((DISK_GB < 10 ? 10 : DISK_GB))
DISK_GB=$((DISK_GB > 200 ? 200 : DISK_GB))
# RAM: 0.5x data + 2GB, minimum 4GB, max 64GB
RAM_GB=$((TAR_SIZE_GB / 2 + 2))
RAM_GB=$((RAM_GB < 4 ? 4 : RAM_GB))
RAM_GB=$((RAM_GB > 64 ? 64 : RAM_GB))
echo "Tar size: ${TAR_SIZE_GB}GB, requesting Disk: ${DISK_GB}GB, RAM: ${RAM_GB}GB"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Submit job and pipe output
output=$( \
  jobsub_submit --disk=${DISK_GB}GB --expected-lifetime=2h --memory=${RAM_GB}GB -G coupp \
    --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE,DEDICATED \
	--tar_file_name dropbox:///${SCRIPT_DIR}/${TARBALL} \
	-N 1 \
	file://${SCRIPT_DIR}/gridjob.sh \
	"${TEMP_DIR}/${RUN_ID}.tar" \
    "${OUT_DIR}" \
)
echo "$output"
# Retrieve Job ID from output, and save to list
JOB_ID=$(echo "$output" | grep -oP '\d+\.\d+@\S+\.fnal\.gov')
echo "$(date '+%Y-%m-%d %H:%M:%S'), ${RUN_ID}, ${JOB_ID}" >> "${LIST_DIR}"
echo "Job ${JOB_ID} successfully submitted at $(date '+%Y-%m-%d %H:%M:%S')"
rm $TARBALL
)

if [[ $? -ne 0 ]]; then
  echo "An error occurred. Quitting."
  rm -f "LAr10ana*.tar"
fi
