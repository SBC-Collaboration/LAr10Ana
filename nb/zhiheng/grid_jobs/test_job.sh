#!/bin/bash
OUTDIR=$1
LOG="${CLUSTER}_${PROCESS}.log"
OUT=output

echo "CONFIGURE ENVIRONMENT" >>${LOG} 2>&1
source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
spack load ifdhc@2.7.1%gcc@11.4.1

echo "INITIALIZE CONDA ENV"
source ${INPUT_TAR_DIR_LOCAL}/setup.sh
export PYTHONPATH=$PYTHONPATH:${INPUT_TAR_DIR_LOCAL}
echo $PYTHONPATH

# echo "COPY RUN TO LOCAL DIRECTORY"
# LOCAL_RUNDIR=`basename $RUNDIR`
# ifdh cp -r $RUNDIR $LOCAL_RUNDIR >>${LOG} 2>&1

python ${INPUT_TAR_DIR_LOCAL}/zhiheng/test_job.py

ifdh mkdir ${OUTDIR}/${CLUSTER}_${PROCESS}
ifdh cp -D ${LOG} ${OUT}/* ${OUTDIR}/${CLUSTER}_${PROCESS}/
ifdh_exit_code=$?

exit ${ifdh_exit_code}
