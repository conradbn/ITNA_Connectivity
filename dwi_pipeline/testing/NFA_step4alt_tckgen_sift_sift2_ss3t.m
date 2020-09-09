%% Run tractography using WM FODs from the SS3T/Dhollander pipeline
% Run track generation algorithm with 10 million streamlines, then
% additonal streamlines from the regions of interest. Combine these
% streamlines and then run SIFT2 streamline weighting process
setenv('PATH', [getenv('PATH') ':/usr/local/bin']); % Required on nbl_imac2 when starting from Applications (not terminal)
purge;
local_dir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData';
cd(start_dir);
sub_dirs = dir('Price_*');

for ii = 1%:numel(sub_dirs)
    cd(start_dir)
    cd([sub_dirs(ii).name '/PREPROCESSED']);
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{2};
%    
%     % Copy data to local directory, then decompress
    local_dir_sub = [local_dir '/' sub];
%     mkdir(local_dir_sub);
%     unix(['cp -f fod_wm_norm.mif.gz ' local_dir_sub]);
%     unix(['cp -f aseg_5ttgen.mif ' local_dir_sub]);
%     unix(['cp -f aseg_5tt_gmwmi.mif ' local_dir_sub]);
%     unix(['cp -f roi_PP19_MNI152/*.binary_vol_al2dwi.nii.gz ' local_dir_sub]);
    cd(local_dir_sub)
%     
%     unix('pigz -d fod_wm_norm.mif.gz');
%     
%     % Generate streamlines for wholebrain
%     % NOTE - the algorithm and cutoff settings are actually the defaults,
%     % but I set them here for clarity
%     unix(['tckgen fod_wm_norm.mif tracks_ss3t_' sub '.tck'...
%       ' -backtrack -crop_at_gmwmi -info'...
%       ' -act aseg_5ttgen.mif '...
%       ' -seed_gmwmi aseg_5tt_gmwmi.mif'...
%       ' -algorithm iFOD2'...
%       ' -select 10M'...
%       ' -cutoff 0.05']);
%     
%     % Generate streamlines for ROIs
%     % NOTE - the algorithm and cutoff settings are actually the defaults,
%     % but I set them here for clarity
    rois = dir('*binary_vol_al2dwi.nii.gz');
%     for jj = 1:numel(rois)
%         r = rois(jj).name;
%         unix(['tckgen fod_wm_norm.mif tracks_ss3t_' sub '_' r(1:8) '.tck'...
%           ' -backtrack -crop_at_gmwmi -info'...
%           ' -act aseg_5ttgen.mif '...
%           ' -seed_gmwmi ' r...
%           ' -algorithm iFOD2'...
%           ' -select 500k'...
%           ' -cutoff 0.05']);
%     end
%     
% %     % Move wholebrain streamlines to subject directory
% %     unix(['mv ../tracks_ss3t_' sub '.tck ' local_dir_sub]);
%     
%     % Combine streamlines 
%     unix(['tckedit *' sub '*.tck tracks_ss3t_' sub '_combined.tck']);

    %% Run original SIFT algorithm to partially filter streamlines based on diffusion signal
    % Stop when 10 million streamlines remain
    unix(['tcksift tracks_ss3t_' sub '_combined.tck'...
         ' fod_wm_norm.mif.gz tracks_ss3t_' sub '_combined_sift10M.tck'...
         ' -act aseg_5ttgen.mif.gz -info '...
         ' -term_number 10M']);
        
    %% Run SIFT2 algorithm to weight streamlines based on diffusion signal, such
    % that they are more closely related to true fiber densities. This allows
    % for more valid quantitative tractography measures.
    % NOTE - to make use of SIFT2 output, the weights txt file must be included
    % in further track quantitification tools using "-tck_weights_in" flag
    % ALSO NOTE - the I removed the -fd_scale_gm flag due to the use of
    % multi-tissue FOD estimation
    unix(['tcksift2 tracks_ss3t_' sub '_combined_sift10M.tck'...
         ' fod_wm_norm.mif.gz tracks_sift2_weights_ss3t_' sub '_combined_sift10M.txt'...
         ' -act aseg_5ttgen.mif.gz -info '...
         ' -out_mu tracks_sift2_weights_ss3t_' sub '_combined_sift10M_prop_coeff.txt']);
    
     
    %% Create the reduced streamline maps for the ROIs 
    
    % Remove headers from weight files (ONLY IF previous implemented by ss3t install)
    %unix(['echo "$(tail -n +2 tracks_sift2_weights_ss3t_' sub '_combined.txt)" > tracks_sift2_weights_ss3t_' sub '_combined.txt']);
    
    for jj = 1:numel(rois)
        r = rois(jj).name;
        % Only look for streamlines with an end in ROI
        unix(['tckedit -force tracks_ss3t_' sub '_combined.tck -ends_only'...
            ' -tck_weights_in tracks_sift2_weights_ss3t_' sub '_combined_sift10M.txt'...
            ' -tck_weights_out tracks_sift2_weights_ss3t_' sub '_combined_sift10M_' r(1:8) '.txt'...
            ' -include ' r(1:8) '.binary_vol_al2dwi.nii.gz '...
            ' tracks_ss3t_' sub '_combined_sift10M_' r(1:8) '.tck']);
    end
    
           
%     % Remove the input files (LEAVE THERE FOR NOW)
%     unix(['rm -f fod_wm_norm.mif aseg_5ttgen.mif aseg_5tt_gmwmi.mif '...
%          ' *.binary_vol_al2dwi.nii.gz']);    
    
end