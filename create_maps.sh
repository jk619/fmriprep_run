#! /bin/bash

## the following script can be used to run analyzePRF on the HPC cluster
## using the command 'sbatch runPrincePRF.sh'

## First we set up some SBATCH directives. Note that these are hard values
## if you go beyond your job will be killed by SLURM

#SBATCH --job-name=create_maps
#SBATCH -a 0  # run this script as 2 jobs with SLURM_ARRAY_TASK_ID = 0 and 1. Add more numbers for more jobs!
#SBATCH --nodes=1 # nodes per job
#SBATCH --cpus-per-task=4 #~2 days to run PRFs
#SBATCH --mem=16g # More memory you request the less priority you get
#SBATCH --time=00:15:00 # Max request to be safe...
#SBATCH --output=/scratch/jk7127/logs/create_maps_out_%x-%a.txt # Define output log location
#SBATCH --error=/scratch/jk7127/logs/create_maps_err_%x-%a.txt # and the error logs for when it inevitably crashes
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
warning off
subject = subject(1:end-1)
disp(subject)

bidsfolder = '/scratch/jk7127/NEI/BIDS/'
fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer', ['sub-' subject]);
prfpath = sprintf('%sderivatives/prfvista/sub-%s/ses-nyu3t01/',bidsfolder,subject);

mgz = MRIread(fullfile(fspth, 'mri', 'orig.mgz'));


lcurv = read_curv(fullfile(fspth, 'surf', 'lh.curv'));
rcurv = read_curv(fullfile(fspth, 'surf', 'rh.curv'));


leftidx  = 1:numel(lcurv);
rightidx = (1:numel(rcurv))+numel(lcurv);

hemi{1} = zeros(length(lcurv),1);
hemi{2} = zeros(length(rcurv),1);

%%


d = dir(sprintf('%s*results*',prfpath));
files = natsort({d.name});

ssigma = [];
vvexpl = [];
aangle = [];
aangle_adj = [];
eeccen = [];
xx = [];
yy = [];

for f = 1 : length(files)

    tmp = load(sprintf('%s/%s',prfpath,files{f}));


  mysigma = tmp.results.model{1}.sigma.major;
    myvexpl = 1 - (tmp.results.model{1}.rss ./ tmp.results.model{1}.rawrss);
    myangle = atan2(-tmp.results.model{1}.y0,tmp.results.model{1}.x0);
    myangle_adj = (mod(90 - 180/pi * myangle + 180, 360) - 180);

    myx     = tmp.results.model{1}.x0;
    myy     = tmp.results.model{1}.y0;
    myeccen =  sqrt(tmp.results.model{1}.x0.^2+tmp.results.model{1}.y0.^2);

    ssigma = [ssigma mysigma];
    vvexpl = [vvexpl myvexpl];
    aangle = [aangle myangle];
    aangle_adj = [aangle_adj myangle_adj];
    xx = [xx myx];
    yy = [yy myy];
    eeccen = [eeccen myeccen];

end

mgz.vol = aangle(leftidx);
MRIwrite(mgz, fullfile(prfpath, 'lh.angle.mgz'));
mgz.vol = aangle(rightidx);
MRIwrite(mgz, fullfile(prfpath, 'rh.angle.mgz'));


mgz.vol = aangle_adj(leftidx);
MRIwrite(mgz, fullfile(prfpath, 'lh.angle_adj.mgz'));

mgz.vol = aangle_adj(rightidx);
MRIwrite(mgz, fullfile(prfpath, 'rh.angle_adj.mgz'));



mgz.vol = eeccen(leftidx);
MRIwrite(mgz, fullfile(prfpath, 'lh.eccen.mgz'));
mgz.vol = eeccen(rightidx);
MRIwrite(mgz, fullfile(prfpath, 'rh.eccen.mgz'));

mgz.vol = ssigma(leftidx);
MRIwrite(mgz, fullfile(prfpath, 'lh.sigma.mgz'));
mgz.vol = ssigma(rightidx);
MRIwrite(mgz, fullfile(prfpath, 'rh.sigma.mgz'));

% r2 (convert from percentage to fraction)
mgz.vol = vvexpl(leftidx);
MRIwrite(mgz, fullfile(prfpath, 'lh.vexpl.mgz'));
mgz.vol = vvexpl(rightidx);
MRIwrite(mgz, fullfile(prfpath, 'rh.vexpl.mgz'));

mgz.vol = xx(leftidx);
MRIwrite(mgz, fullfile(prfpath, 'lh.x.mgz'));
mgz.vol = xx(rightidx);
MRIwrite(mgz, fullfile(prfpath, 'rh.x.mgz'));

mgz.vol = yy(leftidx);
MRIwrite(mgz, fullfile(prfpath, 'lh.y.mgz'));
mgz.vol = yy(rightidx);
MRIwrite(mgz, fullfile(prfpath, 'rh.y.mgz'));


EOF

exit 0

