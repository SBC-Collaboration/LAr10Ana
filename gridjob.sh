RUNDIR=$1
OUTDIR=$2
LOG="${CLUSTER}_${PROCESS}.log"
# LOG=/dev/stdout
OUT=output

echo "CONFIGURE ENVIRONMENT" >>${LOG} 2>&1
source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
spack load ifdhc@2.7.1%gcc@11.4.1

echo "INITIALIZE CONDA ENV"
source ${INPUT_TAR_DIR_LOCAL}/setup.sh

export PYTHONPATH=$PYTHONPATH:${INPUT_TAR_DIR_LOCAL}
echo $PYTHONPATH

echo "COPY RUN TO LOCAL DIRECTORY"
LOCAL_RUNDIR=`basename $RUNDIR`
ifdh cp -r $RUNDIR $LOCAL_RUNDIR >>${LOG} 2>&1
extension="${LOCAL_RUNDIR##*.}"

if [[ "$extension" == ".tar" ]]; then
    echo "Extracting .tar file" >>${LOG} 2>&1
    tar -xf $LOCAL_RUNDIR
    LOCAL_RUNDIR=extension
    echo "Processing $LOCAL_RUNDIR" >>${LOG} 2>&1
elif [[ "$extension" == ".tar.gz" ]]; then
    echo "Extracting .tar.gz file" >>${LOG} 2>&1
    tar -xvf $LOCAL_RUNDIR
    LOCAL_RUNDIR=extension
    echo "Processing $LOCAL_RUNDIR" >>${LOG} 2>&1
fi

echo "RUN EVENT DEALER" >>${LOG} 2>&1
python3 ${INPUT_TAR_DIR_LOCAL}/RunEventDealer.py ./$LOCAL_RUNDIR ./${OUT} >>${LOG} 2>&1
echo "EVENT DELAER COMPLETED WITH EXIT CODE $?" >>${LOG} 2>&1 

ifdh mkdir ${OUTDIR}/${CLUSTER}_${PROCESS}
ifdh cp -D ${LOG} ${OUT}/* ${OUTDIR}/${CLUSTER}_${PROCESS}/
ifdh_exit_code=$?

exit ${ifdh_exit_code}
