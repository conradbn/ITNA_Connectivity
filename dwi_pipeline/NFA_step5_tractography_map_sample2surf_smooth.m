%% Postprocessing of tractography
% Create TDI endpoint maps in subject native DWI-space
% Sample maps to subject's standard surface
% Smooth maps on the surface


purge

% NOT HAVING TO DO THIS ON MY MACBOOK PRO FOR SOME REASON (MATLAB STARTED
% FROM TERMINAL)
% % Add MRtrix3 binaries to path 
% setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
% % Add MRtrix3 binaries to path 
% add_AFNI_to_path
% % Dynamic library setup required for running MRtrix and AFNI in same script
% %setenv('DYLD_LIBRARY_PATH',[ '/usr/lib:/usr/local/lib:' getenv('DYLD_LIBRARY_PATH')]);
% setenv('DYLD_FALLBACK_LIBRARY_PATH',[ '~/abin:' getenv('DYLD_FALLBACK_LIBRARY_PATH') ])

%% Setup
% Set directories
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData';
cd(start_dir);
sub_dirs = dir('Price_*');
local_dir = '/Volumes/BensHD_2020//Price_NFA_Tractography_MNI152';

% % Transforms
affine_dwi = 'brain.finalsurfs_al2dwi_mat.aff12.1D';

%% Loop
for ii = 1:numel(sub_dirs)
    cd(start_dir)
    cd([sub_dirs(ii).name '/PREPROCESSED/warp_T1_to_MNI']);
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{2};
    disp(['***** WORKING ON ' sub ' ******'])
    % Copy to temp directory
    tmp_dir = [local_dir '/' sub];
    cd(tmp_dir);
    
    % SUMA dir
    suma_dir = ['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sub '_proc/' sub '.freesurfer/SUMA'];
    
    %% Create TDI maps in subject native space
    rois = dir('*binary_vol_al2dwi.nii.gz');
    for jj = 1:numel(rois)
        r = rois(jj).name(1:8);
        if strcmp(r(1),'.') % Skip files starting with .
            continue
        else
            % Map track file to endpoint density image
            tck = ['tracks_ss3t_' sub '_combined_' r];
            unix(['tckmap -force -template brain.finalsurfs_al2dwi.nii.gz -ends_only -backtrack -tck_weights_in'...
                ' tracks_sift2_weights_ss3t_' sub '_combined_' r '.txt'...
                ' ' tck '.tck'...
                ' ' tck '.TDI_ends.nii']);
            % Scale by the subject-specific proportionality coefficient
            mu = load(['tracks_sift2_weights_ss3t_' sub '_combined_prop_coeff.txt']);
            unix(['3dcalc -a ' tck '.TDI_ends.nii'...
                ' -prefix ' tck '.TDI_ends.norm.nii'...
                ' -expr "a* ' num2str(mu) '"']);
            % Transform TDI endpoint map to anatomical space
            unix(['cat_matvec ' affine_dwi ' -I > inv_' affine_dwi]);
            unix(['3dAllineate -1Dmatrix_apply inv_' affine_dwi...
                ' -master ' suma_dir '/' sub '_SurfVol.nii'...
                ' -input ' tck '.TDI_ends.norm.nii'...
                ' -prefix ' tck '.TDI_ends.norm.al2anat.nii']);
            % For each hemisphere
            hemi = {'lh','rh'};
            for hh = 1:2
                h = hemi{hh};
                % Sample the TDI endpoint map to the surface
                unix(['3dVol2Surf -overwrite -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
                    ' -sv ' suma_dir '/' sub '_SurfVol.nii'...
                    ' -surf_A smoothwm'...
                    ' -f_index nodes'...
                    ' -map_func mask'...
                    ' -oob_value 0'...
                    ' -grid_parent ' tck '.TDI_ends.norm.al2anat.nii'...
                    ' -out_niml ' tck '.TDI_ends.norm.al2anat.' h '.niml.dset']);
                    % ' -use_norms -norm_length -1'...
                    % ' -f_steps 10'...
                % Smooth on surface and log transform data
                k = '6';
                unix(['SurfSmooth -overwrite -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
                    ' -surf_A smoothwm -met HEAT_07 -fwhm ' k...
                    ' -input ' tck '.TDI_ends.norm.al2anat.' h '.niml.dset'...
                    ' -output ' tck '.TDI_ends.norm.al2anat.' h '.' k 'mm.niml.dset']);
                unix(['3dcalc -overwrite -prefix ' tck '.TDI_ends.norm.al2anat.' h '.' k 'mm.log.niml.dset'...
                     ' -a ' tck '.TDI_ends.norm.al2anat.' h '.' k 'mm.niml.dset'...
                     ' -expr "log(a)"'])
                k = '4';
                unix(['SurfSmooth -overwrite -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
                    ' -surf_A smoothwm -met HEAT_07 -fwhm ' k...
                    ' -input ' tck '.TDI_ends.norm.al2anat.' h '.niml.dset'...
                    ' -output ' tck '.TDI_ends.norm.al2anat.' h '.' k 'mm.niml.dset']);
                unix(['3dcalc -overwrite -prefix ' tck '.TDI_ends.norm.al2anat.' h '.' k 'mm.log.niml.dset'...
                     ' -a ' tck '.TDI_ends.norm.al2anat.' h '.' k 'mm.niml.dset'...
                     ' -expr "log(a)"'])
            end
               
        end
    end
    %Compress all the new nifti files
    unix('pigz *.nii');
end