#! /bin/bash

## the following script can be used to run analyzePRF on the HPC cluster
## using the command 'sbatch runPrincePRF.sh'

## First we set up some SBATCH directives. Note that these are hard values
## if you go beyond your job will be killed by SLURM

#SBATCH --job-name=create-dataset
#SBATCH -a 0  # run this script as 2 jobs with SLURM_ARRAY_TASK_ID = 0 and 1. Add more numbers for more jobs!
#SBATCH --nodes=1 # nodes per job
#SBATCH --cpus-per-task=4 #~2 days to run PRFs
#SBATCH --mem=16g # More memory you request the less priority you get
#SBATCH --time=00:10:00 # Max request to be safe...
#SBATCH --output=/scratch/jk7127/logs/create_data_out_%x-%a.txt # Define output log location
#SBATCH --error=/scratch/jk7127/logs/create_data_err_%x-%a.txt # and the error logs for when it inevitably crashes
#SBATCH --mail-user=jk7127@nyu.edu #email
#SBATCH --mail-type=END #email me when it crashes or better, ends



if [[ $# -eq 0 ]] ; then
    echo 'You need to provice Subject ID'
    exit 0
fi


# load matlab module
module load matlab/2020b
module load freesurfer/6.0.0
source setup.sh $1


matlab -nodesktop -nodisplay -nosplash <<EOF
[~,subject] = system('echo $SUBJECT_ID')
tbUse vistasoft
subject = subject(1:end-1)
disp(subject)

hemi = {'L';'R'}
mypath = sprintf('/scratch/jk7127/NEI/BIDS/derivatives/fmriprep/sub-%s/ses-nyu3t01/func',subject)

tasks = {'bar1';'bar2';'bar3';'wedgering1';'wedgering2';'wedgering3'}

for t = 1 : length(tasks)
for h = 1 : length(hemi)

    d = dir(sprintf('%s/*task-%s_space-fsnative_hemi-%s_bold.func.gii',mypath,tasks{t},hemi{h}))


        [~, fname] = fileparts(d.name);
       mynames{t,h} = sprintf('%s/%s.mgz',mypath,fname)
        str = sprintf('mri_convert %s/%s.gii %s/%s.mgz',mypath,fname, mypath,fname);
        system(str);


end
end


datafiles_all = cell(length(mynames),1);

for f = 1 : length(mynames)
    
    data = [];
    for h = 1 : length(hemi)
        
        tmp = MRIread(mynames{f,h});
        data = cat(1,data,squeeze(tmp.vol));
        
    end
    datafiles_all{f} = data;
end
prfpath = sprintf('/scratch/jk7127/NEI/BIDS/derivatives/prfvista/sub-%s/ses-nyu3t01/',subject)

mkdir(prfpath)

datafiles{1} = nanmean(cat(3,datafiles_all{1},datafiles_all{2},datafiles_all{3}),3);
datafiles{2} = nanmean(cat(3,datafiles_all{4},datafiles_all{5},datafiles_all{6}),3);
save(sprintf('%s/datafiles.mat',prfpath),'datafiles','-v7.3')


warning off


EOF

exit 0
