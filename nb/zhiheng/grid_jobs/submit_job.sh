#!/bin/bash

LAR10ANA_DIR="/exp/e961/app/users/zsheng/LAr10Ana/"
CURRENT_DIR=$(pwd)
shopt -s nullglob

cd ${LAR10ANA_DIR}
mkdir -p zhiheng
cp ${CURRENT_DIR}/{*.py,*.sh} zhiheng/
tar -cvf ${CURRENT_DIR}/lar10ana.tar --exclude='*.pyc' \
  ana zhiheng *.py *.sh
rm -r zhiheng/
cd ${CURRENT_DIR}
echo "Successfully created LAr10Ana tar file"

jobsub_submit --disk=1GB --expected-lifetime=10m --memory=1GB -G coupp --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE,DEDICATED \
	--tar_file_name dropbox:///${CURRENT_DIR}/lar10ana.tar \
	-N 1 \
	file://${CURRENT_DIR}/test_job.sh \
	/pnfs/coupp/scratch/zhiheng/