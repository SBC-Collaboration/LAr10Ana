tar -cf LAr10ana.tar --exclude='*.pyc'  *.py ana setup.sh
echo "LAr10ana is tarred"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
RUN_ID="20251119_9"
DATA_PATH="/pnfs/coupp/storage/SBC-25-data/${RUN_ID}.tar"
OUTPUT_DIR="/pnfs/coupp/scratch/zhiheng/grid_output"

output=$( \
jobsub_submit --disk=50GB --expected-lifetime=8h --memory=6GB -G coupp \
    --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE,DEDICATED \
	--tar_file_name dropbox:///${SCRIPT_DIR}/LAr10ana.tar \
	-N 1 \
	file://${SCRIPT_DIR}/gridjob.sh \
	${DATA_PATH} ${OUTPUT_DIR} \
)
echo "$output"
JOB_ID=$(echo "$output" | grep -oP '\d+\.\d+@\S+\.fnal\.gov')
echo "Job ${JOB_ID} is submitted."
echo "$(date '+%Y-%m-%d %H:%M:%S'), ${RUN_ID}, ${JOB_ID}" >> "${OUTPUT_DIR}/logs/job_list.csv"
rm LAr10ana.tar