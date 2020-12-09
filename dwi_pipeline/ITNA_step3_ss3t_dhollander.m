%% SS3T CSD
% Run single-shell 3 tissue constrained spherical deconvolution
% Note this currently (5/20/20) requires the MRTrix3 fork called
% MRtrix3Tissue from Dhollander et al. To avoid conflict with MRtrix3
% install on the lab iMacs, I installed/ran this package on my MacBook Pro.

purge;
local_dir = '/Users/benconrad/Desktop/dwi_temp';
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData';
cd(start_dir);
sub_dirs = dir('Price_*');

% Set group level tissue response functions
wm_rsp = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData/dhollander_wm_response_group_mean.txt';
gm_rsp = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData/dhollander_gm_response_group_mean.txt';
csf_rsp = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData/dhollander_csf_response_group_mean.txt';

for ii = 1:numel(sub_dirs)
    cd(start_dir)
    cd([sub_dirs(ii).name '/PREPROCESSED']);
    
    % Compress/Copy data to local directory, then decompress
    unix(['cp -f dwmri_pre_biascorr.mif ' local_dir]);
    cd(local_dir)
%     unix('pigz -k dwmri_pre_biascorr.mif');
%     unix(['mv -f dwmri_pre_biascorr.mif.gz ' local_dir]);
%     cd(local_dir)
%     unix('pigz -d dwmri_pre_biascorr.mif.gz');
    
    % Get dwi mask
    unix('dwi2mask dwmri_pre_biascorr.mif dwmri_mask.mif')
    
    % Run single-shell 3 tissue constrained spherical deconvolution to get
    % FODs (ss3t_csd)
    unix(['/Users/benconrad/Downloads/MRtrix3Tissue-3Tissue_v5.2.8/bin/ss3t_csd_beta1 dwmri_pre_biascorr.mif '...
          wm_rsp ' fod_wm.mif '...
          gm_rsp ' fod_gm.mif '...
          csf_rsp ' fod_csf.mif '...
          ' -mask dwmri_mask.mif']);

    % Run the normalizaton/bias correction step
    unix(['/Users/benconrad/Downloads/MRtrix3Tissue-3Tissue_v5.2.8/bin/mtnormalise '...
         ' fod_wm.mif fod_wm_norm.mif'...
         ' fod_gm.mif fod_gm_norm.mif'...
         ' fod_csf.mif fod_csf_norm.mif'...
         ' -mask dwmri_mask.mif']); 
    
    % Send data back to the server
    % First remove the raw dwi data since it is already on server
    unix('rm -f dwmri_pre_biascorr.mif');
    unix('pigz *.mif');
    unix(['mv -f ' local_dir '/*.gz ' start_dir '/' sub_dirs(ii).name '/PREPROCESSED/']);
    %unix(['gunzip ' start_dir '/' sub_dirs(ii).name '/PREPROCESSED/*.mif.gz']);
    %unix('rm -f *.mif.gz');
    
end