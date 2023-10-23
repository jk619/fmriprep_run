#! /bin/bash

## the following script can be used to run analyzePRF on the HPC cluster
## using the command 'sbatch runPrincePRF.sh'

## First we set up some SBATCH directives. Note that these are hard values
## if you go beyond your job will be killed by SLURM

#SBATCH --job-name=docker-vista
#SBATCH -a 0  # run this script as 2 jobs with SLURM_ARRAY_TASK_ID = 0 and 1. Add more numbers for more jobs!
#SBATCH --nodes=1 # nodes per job
#SBATCH --cpus-per-task=16 #~2 days to run PRFs
#SBATCH --mem=16g # More memory you request the less priority you get
#SBATCH --time=48:00:00 # Max request to be safe...
#SBATCH --output=/scratch/jk7127/logs/prf_vista_out_%x-%a.txt # Define output log location
#SBATCH --error=/scratch/jk7127/logs/prf_vista_err_%x-%a.txt # and the error logs for when it inevitably crashes
#SBATCH --mail-user=jk7127@nyu.edu #email
#SBATCH --mail-type=END #email me when it crashes or better, ends

# load matlab module
module load matlab/2020b
module load freesurfer/6.0.0
source setup.sh $1


matlab -nodesktop -nodisplay -nosplash <<EOF
[~,subject] = system('echo $SUBJECT_ID')
tbUse vistasoft
subject = subject(1:end-1)
disp(subject)


mypath = sprintf('/scratch/jk7127/NEI/BIDS/derivatives/fmriprep/sub-%s/ses-nyu3t01/func',subject)


hemi = {'L';'R'}
% convert to mgz using freesurfer
%for h = 1 : length(hemi)
%    
%    d = dir(sprintf('%s/*fsnative_hemi-%s_bold.func.gii',mypath,hemi{h}))
%    
%    for ii = 1:length(d)
%        [~, fname] = fileparts(d(ii).name);
%        mynames{ii,h} = sprintf('%s/%s.mgz',mypath,fname)
%        str = sprintf('mri_convert %s/%s.gii %s/%s.mgz',mypath,fname, mypath,fname);
%        system(str);
%    end
%    
%end
%%
%datafiles = cell(length(mynames),1);

%for f = 1 : length(mynames)
%    
%    data = [];
%    for h = 1 : length(hemi)
%        
%        tmp = MRIread(mynames{ii,h});
%        data = cat(1,data,squeeze(tmp.vol));
%        
%    end
%    datafiles{f} = data;
%end
prfpath = sprintf('/scratch/jk7127/NEI/BIDS/derivatives/prfvista/sub-%s/ses-nyu3t01/',subject)

mkdir(prfpath)

%save(sprintf('%s/datafiles.mat',prfpath),'datafiles','-v7.3')


warning off

% Different subjects have different run numbers so I tend to process them in groups
% that have the same EPI count.

% run the wrapper stript
setenv('MCR_CACHE_ROOT','/scratch/jk7127/cache/')
load stimfiles
run_prfVista();

disp('Complete!');


EOF

exit 0
