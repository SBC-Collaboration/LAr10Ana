tar -cvf LAr10ana.tar *.py ana setup.sh
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
jobsub_submit --disk=50GB --expected-lifetime=8h --memory=6GB -G coupp --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE,DEDICATED \
	--tar_file_name dropbox:///${SCRIPT_DIR}/LAr10ana.tar \
	-N 1 \
	file://${SCRIPT_DIR}/gridjob.sh \
	/pnfs/coupp/scratch/users/gputnam/SBC-25-daqdata/20250611_1 \
	/pnfs/coupp/scratch/users/gputnam/test-acoustic-grid
