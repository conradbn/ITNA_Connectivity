%% Run tractography using WM FODs from the SS3T/Dhollander pipeline
% Run track generation algorithm with 50 million streamlines. Then run
% SIFT2 streamline weighting process. 
%
% NOTES: I found that the distribution of weights was more appropriate in
% this method compared to that when I oversampled (500K) streamlines from
% the ROIs and combined with 10M wholebrain. Nearly all the additional
% streamlines in the oversampling case were downweighted towards 0, and
% could have unintended effects on the regularization procedure performed
% by SIFT2. However, at only 10M streamlines wholebrain, the resulting
% tracking from the ROIs was a bit sparse, and bumping to 50M got this
% closer to the oversampling case. An added benefit should also be that
% there is now increased sampling from the entire cortex, so there is less
% bias in the directionality of connections to the ROIs (as compared to
% explicitly seeding from the ROIs outward). I did notice as well that the
% endpoint density values were inflated in the oversampling case, which is
% worrying and made me further conclude that that approach was introducing
% issues in the weight assignment by SIFT2.

%setenv('PATH', [getenv('PATH') ':/usr/local/bin']); % Required on nbl_imac2 when starting from Applications (not terminal)
purge;
local_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface';
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData';
cd(start_dir);
sub_dirs = dir('Price_*');

for ii = 1:numel(sub_dirs)
    cd(start_dir)
    cd([sub_dirs(ii).name '/PREPROCESSED']);
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{2};
   
    % Copy data to local directory, then decompress
    local_dir_sub = [local_dir '/' sub];
    mkdir(local_dir_sub);
    unix(['cp -f fod_wm_norm.mif.gz ' local_dir_sub]);
    unix(['cp -f aseg_5ttgen.mif ' local_dir_sub]);
    unix(['cp -f aseg_5tt_gmwmi.mif ' local_dir_sub]);
    unix(['cp -f roi_PP19_MNI152/*.binary_vol_al2dwi.nii.gz ' local_dir_sub]);
    cd(local_dir_sub)
    
    unix('pigz -d fod_wm_norm.mif.gz');
    
    % Generate streamlines for wholebrain
    % NOTE - the algorithm and cutoff settings are actually the defaults,
    % but I set them here for clarity
    unix(['tckgen fod_wm_norm.mif tracks_ss3t_' sub '_50M.tck'...
      ' -backtrack -crop_at_gmwmi -info'...
      ' -act aseg_5ttgen.mif '...
      ' -seed_gmwmi aseg_5tt_gmwmi.mif'...
      ' -algorithm iFOD2'...
      ' -select 50M'...
      ' -cutoff 0.05']);
        
    %% Run SIFT2 algorithm to weight streamlines based on diffusion signal, such
    % that they are more closely related to true fiber densities. This allows
    % for more valid quantitative tractography measures.
    % NOTE - to make use of SIFT2 output, the weights txt file must be included
    % in further track quantitification tools using "-tck_weights_in" flag
    % ALSO NOTE - the I removed the -fd_scale_gm flag due to the use of
    % multi-tissue FOD estimation
    unix(['tcksift2 tracks_ss3t_' sub '_50M.tck'...
         ' fod_wm_norm.mif tracks_sift2_weights_ss3t_' sub '_50M.txt'...
         ' -act aseg_5ttgen.mif -info '...
         ' -out_mu tracks_sift2_weights_ss3t_' sub '_50M_prop_coeff.txt']);    
    
end