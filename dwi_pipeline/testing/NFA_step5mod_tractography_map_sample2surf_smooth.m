%% Postprocessing of tractography results
% General outline of processes:
% - Get tracks intersecting the ROIs
% - Create track-weighted maps in subject native DWI-space
% - Transform to anatomical space
% - Sample maps to subject's standard surface
% - Smooth maps on the surface

%% Add binaries to path
% NOT HAVING TO DO THIS ON MY MACBOOK PRO FOR SOME REASON (MATLAB STARTED
% FROM TERMINAL)
% % Add MRtrix3 binaries to path 
% setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
% % Add MRtrix3 binaries to path 
% add_AFNI_to_path
% % Dynamic library setup required for running MRtrix and AFNI in same script
% %setenv('DYLD_LIBRARY_PATH',[ '/usr/lib:/usr/local/lib:' getenv('DYLD_LIBRARY_PATH')]);
% setenv('DYLD_FALLBACK_LIBRARY_PATH',[ '~/abin:' getenv('DYLD_FALLBACK_LIBRARY_PATH') ])

%% Setup to run for all subjects
purge
% Set directories
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData';
cd(start_dir);
sub_dirs = dir('Price_*');
local_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface';

%% Loop through subjects
for ii = 2:numel(sub_dirs)
    
    cd(start_dir)
    
    %% Subject-specific setup
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{2};
    disp(['***** WORKING ON ' sub ' ******'])
   
    % Set the local subject directory
    tmp_dir = [local_dir '/' sub];
%     
%     % Copy anatomical and anat2dwi transform from server
%     unix(['cp ' sub_dirs(ii).name '/PREPROCESSED/brain.finalsurfs_al2dwi.nii.gz ' tmp_dir]);
%     unix(['cp ' sub_dirs(ii).name '/PREPROCESSED/brain.finalsurfs_al2dwi_mat.aff12.1D ' tmp_dir]);
%     
    cd(tmp_dir);
    
    % Set some variables that will be used multiple times
    suma_dir = ['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sub '_proc/' sub '.freesurfer/SUMA'];
    surfvol = [suma_dir '/' sub '_SurfVol.nii'];
    affine_dwi = 'brain.finalsurfs_al2dwi_mat.aff12.1D';
    mu = load(['tracks_sift2_weights_ss3t_' sub '_50M_prop_coeff.txt']);
    
    %% Get streamlines that intersect each ROI
    rois = dir('*binary_vol_al2dwi.nii.gz');
    tck = ['tracks_ss3t_' sub '_50M.tck'];
    tck_weights = ['tracks_sift2_weights_ss3t_' sub '_50M.txt'];
%     for jj = 1:numel(rois)
%         r = rois(jj).name(1:8);
%         if strcmp(r(1),'.') % Skip files starting with .
%             continue
%         else
%             % Only look for streamlines with an end in ROI
%             unix(['tckedit -force -ends_only ' tck...
%                 ' -tck_weights_in ' tck_weights...
%                 ' -tck_weights_out tracks_sift2_weights_ss3t_' sub '_50M_' r '.txt'...
%                 ' -include ' r(1:8) '.binary_vol_al2dwi.nii.gz '...
%                 ' tracks_ss3t_' sub '_50M_' r(1:8) '.tck']);
%         end
%     end
%     
    %% Create TDI (endpoint density) maps in subject native space for each ROI
    for jj = 1:numel(rois)
        r = rois(jj).name(1:8);
        if strcmp(r(1),'.') % Skip files starting with .
            continue
        else
%             % Map track file to endpoint density image
            tck = ['tracks_ss3t_' sub '_50M_' r];
            tck_weights = ['tracks_sift2_weights_ss3t_' sub '_50M_' r '.txt'];
%             unix(['tckmap -force -template brain.finalsurfs_al2dwi.nii.gz -ends_only -backtrack'...
%                 ' -tck_weights_in ' tck_weights...
%                 ' ' tck '.tck'...
%                 ' ' tck '.TDI_ends.nii']);
%             
%             % Scale by the subject-specific proportionality coefficient (mu)
%             unix(['3dcalc -a ' tck '.TDI_ends.nii'...
%                 ' -prefix ' tck '.TDI_ends.norm.nii'...
%                 ' -expr "a* ' num2str(mu) '"']);
%             
%             % Transform to anatomical space (anat2dwi transform will be inverted)
            in =  [tck '.TDI_ends.norm.nii.gz'];
            out = [tck '.TDI_ends.norm.al2anat.nii.gz'];
%             transform_to_anat(affine_dwi,surfvol,in,out);

            % For each hemisphere
            hemi = {'lh','rh'};
            for hh = 1:2
                h = hemi{hh};
                spec = [suma_dir '/std.141.' sub '_' h '.spec'];
                parent = [tck '.TDI_ends.norm.al2anat.nii'];
                out =    [tck '.TDI_ends.norm.al2anat.' h '.niml.dset'];
                log_xform = 'yes';
                sample2surf_and_smooth(surfvol,spec,parent,out,log_xform)
            end
        end
    end
    
    %% Create wholebrain endpoint density map
%     tck = ['tracks_ss3t_' sub '_50M'];
%     tck_weights = ['tracks_sift2_weights_ss3t_' sub '_50M.txt'];
%     unix(['tckmap -force -ends_only -backtrack'...
%         ' -template brain.finalsurfs_al2dwi.nii.gz'...
%         '  -tck_weights_in ' tck_weights...
%         ' ' tck '.tck'...
%         ' ' tck '.wholebrain_TDI_ends.nii']);
%     
%     % Scale by the subject-specific proportionality coefficient (mu)
%     unix(['3dcalc -a ' tck '.wholebrain_TDI_ends.nii'...
%         ' -prefix ' tck '.wholebrain_TDI_ends.norm.nii'...
%         ' -expr "a* ' num2str(mu) '"']);
%     
%     % Transform to anatomical space (anat2dwi transform will be inverted)
%     in =  [tck '.wholebrain_TDI_ends.norm.nii'];
%     out = [tck '.wholebrain_TDI_ends.norm.al2anat.nii'];
%     transform_to_anat(affine_dwi,surfvol,in,out);
%     
%     % Transform to surface for each hemisphere smooth
%     hemi = {'lh','rh'};
%     for hh = 1:2
%         h = hemi{hh};
%         spec = [suma_dir '/std.141.' sub '_' h '.spec'];
%         parent = [tck '.wholebrain_TDI_ends.norm.al2anat.nii'];
%         out =    [tck '.wholebrain_TDI_ends.norm.al2anat.' h '.niml.dset'];
%         log_xform = 'yes';
%         sample2surf_and_smooth(surfvol,spec,parent,out,log_xform)
%     end
% 
%     %% Create wholebrain streamline length maps (at endpoint voxels only)
%     tck = ['tracks_ss3t_' sub '_50M'];
%     tck_weights = ['tracks_sift2_weights_ss3t_' sub '_50M.txt'];
%     unix(['tckmap -force -contrast length -stat_vox mean'...
%         ' -ends_only -backtrack'...
%         ' -template brain.finalsurfs_al2dwi.nii.gz'...
%         '  -tck_weights_in ' tck_weights...
%         ' ' tck '.tck'...
%         ' ' tck '.wholebrain_length_map.nii']);
%     
%     % Transform to anatomical space (anat2dwi transform will be inverted)
%     in =  [tck '.wholebrain_length_map.nii'];
%     out = [tck '.wholebrain_length_map.al2anat.nii'];
%     transform_to_anat(affine_dwi,surfvol,in,out);
%     
%     % Transform to surface for each hemisphere smooth
%     for hh = 1:2
%         h = hemi{hh};
%         spec = [suma_dir '/std.141.' sub '_' h '.spec'];
%         parent = [tck '.wholebrain_length_map.al2anat.nii'];
%         out =    [tck '.wholebrain_length_map.al2anat.' h '.niml.dset'];
%         log_xform = 'no';
%         sample2surf_and_smooth(surfvol,spec,parent,out,log_xform)
%     end
%     
%     %% Create wholebrain target diversity maps
%     % NO SOLUTION FOR THIS IDEA YET!
%     % It may be positively correlated with average length of streamlines
%     % BUT, it is still distinct information, e.g. long streamlines but all
%     % going to same place... Problem is how to specify "targets" 1.
%     % FreeSurfer regions (but then biased by region size). 2. Same sized
%     % ROIs, e.g. Craddock parcellations. 3. Average distance of
%     % streamline endpoints. But then have to define "distance".
%     % IN ANY CASE, this is a tough algorithmic challenge and would likely
%     % be very computationally costly.
%     
%     %% Create ROI to FreeSurfer parcellation connectivity vectors
%     % This is in order to facilitate between-hemisphere comparisons
%     % Requires AFNI's updated FreeSurfer parcellation renumbering command
%     % (@SUMA_renumber_FS)
%     
%     % Copy parcellation files we need
%     unix(['cp ' suma_dir '/aparc.a2009s+aseg_REN_gmrois.nii.gz ' tmp_dir]); 
%     unix(['cp ' suma_dir '/aparc.a2009s+aseg_REN_all.niml.lt ' tmp_dir]);
% 
%     % Transform parcellation to DWI space
%     unix(['3dAllineate -interp NN -1Dmatrix_apply ' affine_dwi...
%         ' -master ' surfvol...
%         ' -input aparc.a2009s+aseg_REN_gmrois.nii.gz'...
%         ' -prefix aparc.a2009s+aseg_REN_gmrois.al2dwi.nii.gz']);
%     
%     % Loop through ROIs
%     for jj = 1:numel(rois)
%         r = rois(jj).name(1:8);
%         if strcmp(r(1),'.') % Skip files starting with .
%             continue
%         else
%             tck = ['tracks_ss3t_' sub '_50M_' r];
%             tck_weights = ['tracks_sift2_weights_ss3t_' sub '_50M_' r '.txt'];
%             % Extract the connectivity weights
%             unix(['tck2connectome -vector ' tck '.tck -tck_weights_in ' tck_weights...
%                 ' aparc.a2009s+aseg_REN_gmrois.al2dwi.nii.gz '...
%                 tck '.fingerprint.csv']);
%             
%             % Scale by the subject-specific proportionality coefficient (mu)
%             f = readmatrix([tck '.fingerprint.csv'],'NumHeaderLines',1);
%             f = f.*mu;
%             writematrix(f,[tck '.fingerprint.norm.csv']);
%         end
%     end
%     
%     %% Compress all the new nifti files
%     unix('pigz -f *.nii');
%     
end



%% GLOBAL FUNCTIONS 
%% Transform track weighted image from DWI to anatomical space
function [] = transform_to_anat(affine_dwi,surfvol,in,out)
    % Transform TDI endpoint map to anatomical space
    unix(['cat_matvec ' affine_dwi ' -I > inv_' affine_dwi]);
    unix(['3dAllineate -1Dmatrix_apply inv_' affine_dwi...
        ' -master ' surfvol...
        ' -input ' in...
        ' -prefix ' out]);
end

%% Sampling to surface, smooth, also create log transformed version
function [] = sample2surf_and_smooth(surfvol,spec,parent,out,log_xform)
    % Sample the TDI endpoint map to the surface
    unix(['3dVol2Surf -overwrite -spec ' spec...
        ' -sv ' surfvol ...
        ' -surf_A smoothwm'...
        ' -surf_B pial'...
        ' -f_p1_mm -1.0'... % Find max (absolute) value on line segment from pial to wm surface + 1mm inward
        ' -f_steps 10'...
        ' -f_index nodes'...
        ' -map_func max_abs'...
        ' -oob_value 0'...
        ' -grid_parent ' parent...
        ' -out_niml /Users/nbl_imac2/Desktop/temp_out/' out]);
    % Smooth on surface
    k = '6';
    unix(['SurfSmooth -overwrite -spec ' spec...
        ' -surf_A smoothwm -met HEAT_07 -target_fwhm ' k... % -Niter 11 -sigma 0.3679 = params for 4mm'...
        ' -input /Users/nbl_imac2/Desktop/temp_out/' out...
        ' -output /Users/nbl_imac2/Desktop/temp_out/' out(1:end-10) '.' k 'mm.niml.dset']);
    
    if strcmp(log_xform,'yes')
        % Log transform data (to be used in some linear regression
        % analyses)
        unix(['3dcalc -overwrite -prefix /Users/nbl_imac2/Desktop/temp_out/' out(1:end-10) '.' k 'mm.log.niml.dset'...
            ' -a /Users/nbl_imac2/Desktop/temp_out/' out(1:end-10) '.' k 'mm.niml.dset'...
            ' -expr "log(a)"']);
    end
end





