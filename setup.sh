unset PYTHONHOME
unset PYTHONPATH
export PATH=$PATH:/cvmfs/coupp.opensciencegrid.org/LAr10Ana/miniforge3/condabin
source /cvmfs/coupp.opensciencegrid.org/LAr10Ana/miniforge3/etc/profile.d/conda.sh
conda activate env

# add current directory to python path
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export PYTHONPATH=$PYTHONPATH:$SCRIPT_DIR
