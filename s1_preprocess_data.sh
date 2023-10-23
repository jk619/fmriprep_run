## PREPROCESS DATA, including:
#   1. conversion to BIDS
#   2. Defacing
#   3. MRIQC
#   4. FMRIPrep

# Exit upon any error
set -euo pipefail


## Change lines for Study/subject/system. Also add lines for stim files, TSV fles, and eye tracking files. Then try running it on our sample data.
DIRN=`dirname $0`
source $DIRN/setup.sh $1
source "$CODE_DIR"/Subroutines/utils.sh


# Go!
# Sample data from subject wlsubj042, acquired on DATE!?!?!?!
#

###   Global variables:   ###

# Study/subject specific #
dcmFolder="$SAMPLE_DATA_DIR/dicoms"
logFolder=${LOG_DIR}/s1

mkdir -p $logFolder

fsLicense=$(fsLicensePath)


# we'll be running the Docker containers as yourself, not as root:
userID=$(id -u):$(id -g)


###   Get conatiner images:   ###
#container_pull nipreps/fmriprep:20.2.3


# Set up some derived variables that we'll use later:
fsLicenseBasename=$(basename $fsLicense)
fsLicenseFolder=${fsLicense%$fsLicenseBasename}

scp  ${SINGULARITY_PULLFOLDER}/fmriprep_23.1.2.sif  ${SINGULARITY_PULLFOLDER}/fmriprep_23.1.2.sif_$1.sif

###   fMRIPrep:   ###
# fmriprep folder contains the reports and results of 'fmriprep'
container_run  \
           $STUDY_DIR:/data \
           ${fsLicenseFolder}:/FSLicenseFolder:ro \
           nipreps/fmriprep:20.2.3 \
           ${SINGULARITY_PULLFOLDER}/fmriprep_23.1.2.sif_$1.sif \
               "/data \
               /data/derivatives \
               participant \
               --fs-license-file /FSLicenseFolder/$fsLicenseBasename \
               --output-space T1w:res-native fsnative fsaverage MNI152NLin6Asym:res-native \
               --participant_label ${SUBJECT_ID} \
               --skip_bids_validation \
               --omp-nthreads 32 \
	       --mem_mb 64000 \
	       --output-layout legacy \
               --no-submm-recon" \
           
rm ${SINGULARITY_PULLFOLDER}/fmriprep_20.2.3_$1.sif
