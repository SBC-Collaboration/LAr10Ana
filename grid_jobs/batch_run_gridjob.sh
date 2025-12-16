#!/bin/bash
(
set -e

DATA_DIR="/exp/e961/data/SBC-25-daqdata"
RECON_DIR="/pnfs/coupp/persistent/"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo "Starting batch submission of grid jobs..."

# Count total tar files
REGEX="${DATA_DIR}"/*.tar
total=$(ls -1 ${REGEX} 2>/dev/null | wc -l)

if [ $total -eq 0 ]; then
    echo "No .tar files found in ${DATA_DIR}"
    exit 1
fi
echo "Found ${total} .tar files to process"
count=0

# Loop through all .tar files in DATA_DIR
for tar_file in ${REGEX}; do
    # Extract run_id from filename and submit job
    run_id=$(basename "$tar_file" .tar)
    count=$((count + 1))
    echo -e "\n=========================================="
    echo "Processing ${count}/${total}: Run ${run_id}"
    "${SCRIPT_DIR}/run_gridjob.sh" "$run_id"
done

echo "Batch submission of ${total} jobs complete."
)

if [[ $? -ne 0 ]]; then
    echo "An error occurred during batch submission. Quitting."
fi