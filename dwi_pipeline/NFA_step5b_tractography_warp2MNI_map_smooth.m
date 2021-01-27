%% Postprocessing of tractography
% Transform streamlines to MNI space using outputs from AFNI auto_warp.py 
% Create TDI and TDI endpoint maps
% Smooth the TDI maps

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
local_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface';%'/Users/benconrad/Desktop/Price_NFA_Tractography';%'/Users/nbl_imac2/Desktop/Price_NFA_Tractography';

% Inputs
input = 'brain.finalsurfs_al2dwi.nii.gz';
template = '~/abin/MNI152_2009_template.nii.gz'; 

% Transforms
warp = 'anat.un.aff.qw_WARP.nii';
affine = 'anat.un.aff.Xat.1D';%'anat.un.aff.maskwarp.Xat.1D'%'anat.un.aff.nii.Xaff12.1D';
affine_dwi = 'brain.finalsurfs_al2dwi_mat.aff12.1D';

%% Loop
for ii = 1:numel(sub_dirs)
    cd(start_dir)
    cd([sub_dirs(ii).name '/PREPROCESSED/warp_T1_to_MNI']);
    out_dir = pwd;
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{2};
    disp(['***** WORKING ON ' sub ' ******'])
    % Copy to temp directory
    tmp_dir = [local_dir '/' sub];
    mkdir(tmp_dir)
    unix(['cp ../' input ' ' warp ' ' affine ' ../' affine_dwi ' ' tmp_dir]);
    cd(tmp_dir);
    
    % Initialize warp (ORIGINALLY RAN THIS WITH DWI SPACE IMAGE INSTEAD OF
    % TEMPLATE, THIS WAS INCORRECT AND LED TO WEIRD CROPPING IN SEVERAL SUBJECTS)
    unix(['warpinit -force ' template ' identity_warp[].nii']); %input

    % Apply the affine + nonlinear transformation to the identity warp
    % Here we are performing the inverse transform as this is what is required
    % for transformation of tck files (i.e. streamlines)
    for kk = 0:2
        unix(['~/abin/3dNwarpApply -iwarp -nwarp "' warp ' ' affine ' INV(' affine_dwi ')"'...
              ' -source  identity_warp' num2str(kk) '.nii'...
              ' -master ' input... 
              ' -prefix inv_mrtrix_warp' num2str(kk) '.nii']);
    end

    %% Fix warp (replaces 0,0,0 with NaN,NaN,NaN to fit with MRtrix convention)
    unix('warpcorrect -force inv_mrtrix_warp[].nii inv_mrtrix_warp_corrected.mif')% -marker -1');

    %% Transform MNI template to DWI space for verification
    unix(['mrtransform -force ' template ' -warp inv_mrtrix_warp_corrected.mif template_warped2dwi.mif']);
    
    %% Transform track file, create TDI, and smooth 
    rois = dir('*binary_vol_al2dwi.nii.gz');
    for jj = 3%1:numel(rois)
        r = strrep(rois(jj).name,'.binary_vol_al2dwi.nii.gz','');
        % Apply nonlinear transform to MNI space
        unix(['tcktransform -force tracks_ss3t_' sub '_50M_' r '.tck'...
            ' inv_mrtrix_warp_corrected.mif.gz '... % added .gz as it was previously zipped
            ' tracks_ss3t_' sub '_50M_' r '.MNI.tck']);
        % Create TDI map
        unix(['tckmap -force -template ' template ' -precise -backtrack -tck_weights_in'...
            ' tracks_sift2_weights_ss3t_' sub '_50M_' r '.txt'...
            ' tracks_ss3t_' sub '_50M_' r '.MNI.tck'...
            ' tracks_ss3t_' sub '_50M_' r '.MNI.TDI.nii']);
        % Create TDI Ends map (i.e. only the ends of streamlines)
        unix(['tckmap -force -template ' template ' -ends_only -backtrack -tck_weights_in'...
            ' tracks_sift2_weights_ss3t_' sub '_50M_' r '.txt'...
            ' tracks_ss3t_' sub '_50M_' r '.MNI.tck'...
            ' tracks_ss3t_' sub '_50M_' r '.MNI.TDI_ends.nii']);
        % Smooth the TDI maps
        kernel = {'6'}; %{'2' '4' '6' '8'};
        for kk = 1:numel(kernel)
            unix(['mrfilter -force -fwhm ' kernel{kk}...
                ' tracks_ss3t_' sub '_50M_' r '.MNI.TDI.nii'...
                ' smooth'...
                ' tracks_ss3t_' sub '_50M_' r '.MNI.TDI.' kernel{kk} 'mm.nii']);
            unix(['mrfilter -force -fwhm ' kernel{kk}...
                ' tracks_ss3t_' sub '_50M_' r '.MNI.TDI_ends.nii'...
                ' smooth'...
                ' tracks_ss3t_' sub '_50M_' r '.MNI.TDI_ends.' kernel{kk} 'mm.nii']);
        end
    end

   %% Compress all the new nifti files
   unix('pigz -f *.nii');
   unix('pigz *.mif'); 
   
   %% Copy new files to server
   %unix(['rsync -a --progress ' tmp_dir '/* ' out_dir]);
   %unix(['cp -vf ' tmp_dir '/* ' out_dir]);

end