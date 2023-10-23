%tbUse vistasoft
addpath(genpath('./prfVista'))
%load stimfiles.mat
load(sprintf('/scratch/jk7127/NEI/BIDS/derivatives/prfvista_par/sub-%s/ses-nyu3t01/datafiles.mat',subject))
copyfile(sprintf('/scratch/jk7127/NEI/BIDS/derivatives/prfvista_par/sub-%s/ses-nyu3t01/datafiles.mat',subject),sprintf('/scratch/jk7127/NEI/BIDS/derivatives/prfvista/sub-%s/ses-nyu3t01/datafiles.mat',subject))
copyfile('stimfiles.mat',sprintf('%s/stimfiles.mat',prfpath))
disp('loaded')
stimradius = 12.4;
tr = 1;

disp(subject)
%parpool('local',40,'IdleTimeout',180);
results = prfVistasoft(stimfiles,datafiles, stimradius,'tr',tr,'wsearch','coarse to fine and hrf');
save(sprintf('%s/results.mat',prfpath),'results','-v7.3');

