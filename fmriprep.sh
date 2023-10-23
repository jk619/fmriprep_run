#! /bin/bash

var1="/scratch/jk7127/logs/fmriprep_err_ses_%x-%a_$1.txt"
var2="/scratch/jk7127/logs/fmriprep_out_ses_%x-%a_$1.txt"

sbatch --error=$var1 --output=$var2 run_fmriprep.sh $1
