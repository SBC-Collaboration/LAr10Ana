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
tar -cf LAr10ana.tar --exclude='*.pyc' *.py ana setup.sh
echo "Data copied over. LAr10ana is tarred. Ready for job submission."

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Submit job and pipe output
output=$( \
  jobsub_submit --disk=50GB --expected-lifetime=2h --memory=12GB -G coupp \
    --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE,DEDICATED \
	--tar_file_name dropbox:///${SCRIPT_DIR}/LAr10ana.tar \
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
rm LAr10ana.tar
)

if [[ $? -ne 0 ]]; then
  echo "An error occurred. Quitting."
  rm -f LAr10ana.tar
fi
