%% Response function estimation using the dhollander algorithm. 
% Outputs will be averaged across subjects, and then the average used for
% all subjects 

purge;
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData';
cd(start_dir);
sub_dirs = dir('Price_*');
for ii = 1:numel(sub_dirs)
    cd([sub_dirs(ii).name '/PREPROCESSED']);
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{2};

    % Convert nifti/bvec/bval data to mif format (MRtrix image format), use the
    % data before bias correction (since this step will be performed by
    % mtnormalize later, directly on the FODs)
    nii_dwi_proc = 'dwmri.nii.gz'; 
    unix(['mrconvert -fslgrad dwmri.bvec dwmri.bval ' nii_dwi_proc ' dwmri_pre_biascorr.mif']);

    % Run response function estimation
    unix(['dwi2response dhollander dwmri_pre_biascorr.mif '...
          ' dhollander_wm_response_' sub '.txt'...
          ' dhollander_gm_response_' sub '.txt'...
          ' dhollander_csf_response_' sub '.txt'...
          ' -voxels dhollander_voxels_' sub '.mif -force']);
    cd(start_dir)
end