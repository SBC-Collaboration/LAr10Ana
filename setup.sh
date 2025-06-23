unset PYTHONHOME
unset PYTHONPATH
export PATH=$PATH:/exp/e961/app/users/gputnam/conda/condabin/
source /exp/e961/app/users/gputnam/conda/etc/profile.d/conda.sh
conda activate env

# add current directory to python path
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export PYTHONPATH=$PYTHONPATH:$SCRIPT_DIR
