%% Transform ROIs to MNI space to check degree of overlap across subjects
% Note since tractography was performed in subject-DWI space and ROIs were
% defined on subjects' ld141 standard surface, the ROIs in MNI space are
% not necessary for the primary analyses.

purge
% Set directories
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData';
cd(start_dir);
sub_dirs = dir('Price_*');
local_dir = '/Users/benconrad/Desktop/Price_NFA_Tractography';%'/Users/nbl_imac2/Desktop/Price_NFA_Tractography';

% Inputs
input = 'brain.finalsurfs_al2dwi.nii.gz';
template = '~/abin/MNI152_2009_template.nii.gz'; 

% Transforms
warp = 'anat.un.aff.qw_WARP.nii';
affine = 'anat.un.aff.Xat.1D';%'anat.un.aff.maskwarp.Xat.1D'%'anat.un.aff.nii.Xaff12.1D';
affine_dwi = 'brain.finalsurfs_al2dwi_mat.aff12.1D';



for ii = 1:numel(sub_dirs)
    cd(start_dir)
    cd([sub_dirs(ii).name '/PREPROCESSED/warp_T1_to_MNI']);
    out_dir = pwd;
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{2};
    disp(['***** WORKING ON ' sub ' ******'])
    % Copy to temp directory
    tmp_dir = ['/Users/benconrad/Desktop/Price_NFA_Tractography/' sub];
    %unix(['cp ../' input ' ' warp ' ' affine ' ../' affine_dwi ' ' tmp_dir]);
    
    % Copy the ROIs from the SUMA directory
    suma_dir = ['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sub '_proc/' sub '.freesurfer/SUMA'];
    %unix(['cp ' suma_dir '/*PP19*.nii.gz ' tmp_dir]);
    
    % Go to tmp_dir and get roi filenames
    cd(tmp_dir)
    rois = dir('*PP19*.nii.gz');
    
    % Combine ROIs into single images
    unix(['3dcalc -overwrite -prefix combined_PP19_rois_MNI152_' sub '.nii.gz'...
        ' -a std.141.lh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.nii.gz[1] '...
        ' -b std.141.lh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam.nii.gz[1] '...
        ' -c std.141.rh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.nii.gz[1] '...
        ' -d std.141.rh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam.nii.gz[1] '...
        ' -expr "a + b + c + d"']);
    
    unix(['3dcalc -overwrite -prefix combined_PP19_rois_MNI_N27_' sub '.nii.gz'...
        ' -a std.141.lh.PP19_Dp-Da.votc.inflated.10mm_diam.nii.gz[1] '...
        ' -b std.141.lh.PP19_Lp-La.votc.inflated.10mm_diam.nii.gz[1] '...
        ' -expr "a + b"']);
    
    % Apply the affine + nonlinear transformation to the ROIs
    %     for jj = 1:numel(rois)
    %         roi = rois(jj).name;
    %         unix(['3dNwarpApply -nwarp "' warp ' ' affine '"'...
    %             ' -source ' roi...
    %             ' -master ' template...
    %             ' -prefix MNI_' roi]);
    %     end
    
    unix(['3dNwarpApply -overwrite -interp NN -nwarp "' warp ' ' affine '"'...
        ' -source combined_PP19_rois_MNI152_' sub '.nii.gz'...
        ' -master ' template...
        ' -prefix MNI_combined_PP19_rois_MNI152_' sub '.nii.gz']);
    unix(['3dNwarpApply -overwrite -interp NN -nwarp "' warp ' ' affine '"'...
        ' -source combined_PP19_rois_MNI_N27_' sub '.nii.gz'...
        ' -master ' template...
        ' -prefix MNI_combined_PP19_rois_MNI_N27_' sub '.nii.gz']);
          
end


%% Combining the MNI-transformed rois across the group
cd /Users/benconrad/Desktop/Price_NFA_Tractography

unix(['3dMean -overwrite -prefix MNI_combined_PP19_rois_MNI_N27_group_mean.nii.gz '...
     ' */MNI_combined_PP19_rois_MNI_N27_*.nii.gz']);
 
unix(['3dMean -overwrite -prefix MNI_combined_PP19_rois_MNI152_group_mean.nii.gz '...
     ' */MNI_combined_PP19_rois_MNI152_*.nii.gz']);


















