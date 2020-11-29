%% Copy fROIs to DWI Directory and run structural connectivity processing
% Also, creates FreeSurfer volumetric ROIs before connectivity estimation
% MRTrix processing is done on local copy of data to speed up the
% read/write operations

purge;
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData');
sub_dirs = dir('*_proc');

% Set filenames and labels
% rois = {'std.141.lh.PollackandPrice19_Dp-Da.MNI152.votc.inflated.14mm_diam.nii.gz'
%         'std.141.lh.PollackandPrice19_Lp-La.MNI152.votc.inflated.14mm_diam.nii.gz'
%         'std.141.rh.PollackandPrice19_Dp-Da.MNI152.votc.inflated.14mm_diam.nii.gz'
%         'std.141.rh.PollackandPrice19_Lp-La.MNI152.votc.inflated.14mm_diam.nii.gz'};
% conds = {'Dp-Da.lh'
%          'Lp-La.lh'
%          'Dp-Da.rh'
%          'Lp-La.rh'};
rois = {'std.141.rh.PP19_Dp-Da_math.MNI152.votc.inflated.14mm_diam.nii.gz'};
conds = {'Dp-Da_math.rh'};


%% Loop through subjects and copy/combine/resample ROIs
for ss = 1:numel(sub_dirs)
    % Get subject ID and go to freesurfer SUMA folder
    s = strsplit(sub_dirs(ss).name,'_'); 
    s = s{1,1};
    % Print subject progress
    disp('*************************************************');
    disp(['*************** WORKING ON ' s ' ***************']);
    disp('*************************************************');
    
    % Go to subject's SUMA directory
    suma_dir = ['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' s '_proc/' s '.freesurfer/SUMA'];
    cd(['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' s '_proc/' s '.freesurfer/SUMA'])
    
    % Create ROI output directory in DWI preprocessing folder
    roi_dir = ['/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData/Price_' s '_dwi_proc/PREPROCESSED/roi_PP19_MNI152'];
    mkdir(roi_dir);
    cd(roi_dir);
         
    for rr = 1:numel(rois)
        % Copy roi file
        unix(['cp ' suma_dir '/' rois{rr} ' ' roi_dir]);
        % Get binary version of PollackandPrice19 ROIs
        % (originally output as a two volume dataset from AFNI surf2vol process)
        unix(['3dcalc -overwrite -a ' rois{rr} '[1] -prefix ' conds{rr} '.binary_vol.nii.gz -expr "step(a-0)"']);
        % Transform to DWI space
        % Apply transform to the ROI volumes
        unix(['3dAllineate -overwrite -final NN -1Dmatrix_apply ../brain.finalsurfs_al2dwi_mat.aff12.1D'...
            ' -master ../aseg_al2dwi.nii.gz'...
            ' -input ' conds{rr} '.binary_vol.nii.gz'...
            ' -prefix ' conds{rr} '.binary_vol_al2dwi.nii.gz']);
    end
end







