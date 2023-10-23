## Export the path to the download data

# Where does code live?
CODE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Source the subroutines
source "$CODE_DIR"/Subroutines/utils.sh

# Define workspace variables

# Path to downloaded data
# 	If someone calls this with an argument, use that as the path for
# 	SAMPLE_DATA_DIR (though we'll add DownloadedData to it)

# else, set SAMPLE_DATA_DIR to a directory in the user's scratch
# directory
SAMPLE_DATA_DIR="/scratch/$(whoami)/fmriprep_run"    
SINGULARITY_PULLFOLDER="/scratch/$(whoami)/Singularity" 
mkdir -p $SINGULARITY_PULLFOLDER 

# Path to BIDS directory
STUDY_DIR="$SAMPLE_DATA_DIR/BIDS"
# create the necessary directories, just to be safe
mkdir -p $SAMPLE_DATA_DIR 
mkdir -p $STUDY_DIR

# BIDS specific variables
SUBJECT_ID=$1
SESSION_ID=nyu3t01 
LOG_DIR=${STUDY_DIR}/derivatives/logs/sub-${SUBJECT_ID}

# Which container software to use 
CONTAINER_SOFTWARE=`which_software`

export SUBJECT_ID
export SESSION_ID
export SAMPLE_DATA_DIR
export STUDY_DIR
export LOG_DIR
export CODE_DIR
export CONTAINER_SOFTWARE
export SINGULARITY_PULLFOLDER

echo "*******DEFINE VARIABLES*******"
echo "SUBJECT_ID: $SUBJECT_ID"
echo "SESSION_ID: $SESSION_ID"
echo "SAMPLE_DATA_DIR: $SAMPLE_DATA_DIR"
echo "STUDY_DIR: $STUDY_DIR"
echo "LOG_DIR: $LOG_DIR"
echo "CODE_DIR: $CODE_DIR"
echo "CONTAINER_SOFTWARE: $CONTAINER_SOFTWARE"
echo "SINGULARITY_PULLFOLDER: ${SINGULARITY_PULLFOLDER-empty}"

# debug
#  echo $CLUSTER 
#  echo $(on_cluster)
# exit 0
