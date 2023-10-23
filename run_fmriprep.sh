#! /bin/bash

## the following script can be used to run analyzePRF on the HPC cluster
## using the command 'sbatch runPrincePRF.sh'

## First we set up some SBATCH directives. Note that these are hard values
## if you go beyond your job will be killed by SLURM

#SBATCH --job-name=fmriprep
#SBATCH -a 0  # run this script as 2 jobs with SLURM_ARRAY_TASK_ID = 0 and 1. Add more numbers for more jobs!
#SBATCH --nodes=1 # nodes per job
#SBATCH --cpus-per-task=48 #~2 days to run PRFs
#SBATCH --mem=64g # More memory you request the less priority you get
#SBATCH --time=10:00:00 # Max request to be safe...
#SBATCH --mail-user=jk7127@nyu.edu #email
#SBATCH --mail-type=END #email me when it crashes or better, ends

sub=$1
all_subjects=(wlsubj144 wlsubj145)


# load the freesurfer module
module load freesurfer/6.0.0
source /share/apps/freesurfer/6.0.0/SetUpFreeSurfer.sh

sh s1_preprocess_data.sh $sub

exit 0
EOT
