%% Nonlinear registration of T1 to MNI152 template
% Transforms used for later registration of tractography to template space

purge;
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/ProcessedData';
cd(start_dir);
sub_dirs = dir('Price_*');
tmp_dir = '/Users/nbl_imac2/Desktop/sf_afni_tmp/NFA/auto_warp';

for ii = 1:numel(sub_dirs)
    cd(start_dir)
    cd([sub_dirs(ii).name '/PREPROCESSED']);
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{2};
    disp(['***** WORKING ON ' sub ' ******'])

    % Make output directory
    out_dir =pwd; %[pwd '/SSwarper'];
    %unix(['mkdir ' out_dir]);
    
    % Copy T1 to local directory and go there
    %unix(['cp /Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sub '_proc/' sub '.freesurfer/SUMA/' sub '_SurfVol.nii ' tmp_dir])
    unix(['cp brain.finalsurfs.nii ' tmp_dir])
    
% * BETTER TO JUST COMBINE THE TRANSFORMS LATER
      % Transform original T1 to DWI space
%     unix(['3dAllineate -1Dmatrix_apply brain.finalsurfs_al2dwi_mat.aff12.1D'...
%           ' -master ' tmp_dir '/' sub '_SurfVol.nii'...
%           ' -input ' tmp_dir '/' sub '_SurfVol.nii'...
%           ' -prefix ' tmp_dir '/' sub '_SurfVol_al2dwi.nii']);

% * SSwarper IS TAKING A REALLY LONG TIME/NOT FINISHING, USING auto_warp INSTEAD *
%     % Run SSwarper to nonlinearly register the (DWI-aligned) T1 to template space
%     cd(tmp_dir)
%     unix(['@SSwarper'...
%           ' -input ' sub '_SurfVol.nii'...
%           ' -base MNI152_2009_template_SSW.nii.gz'...
%           ' -subid ' sub...
%           ' -giant_move']);  

    % Run auto_warp.py to nonlinearly register the (DWI-aligned) T1 to template space
    cd(tmp_dir)
    unix(['auto_warp.py'...
          ' -input brain.finalsurfs.nii'...
          ' -base MNI152_2009_template.nii.gz'... 
          ' -skull_strip_input no'...
          ' -unifize_input yes'...
          ' -output_dir warp_T1_to_MNI'...
          ' -qw_opts -lite']);  

    unix('rm -f brain.finalsurfs.nii');
    unix(['cp ~/abin/MNI152_2009_template.nii.gz ' tmp_dir '/warp_T1_to_MNI'])
    
    % Make QC images
    unix(['@snapshot_volreg3 ~/abin/MNI152_2009_template.nii.gz '...
           ' warp_T1_to_MNI/brain.finalsurfs.aw.nii  warp_T1_to_MNI/QC_warp_T1_over_MNI_' sub]);
    unix(['@snapshot_volreg3 warp_T1_to_MNI/brain.finalsurfs.aw.nii '...
           '~/abin/MNI152_2009_template.nii.gz warp_T1_to_MNI/QC_warp_MNI_over_T1_' sub]);
    
    % Copy data to the server and remove 
    unix(['rsync -a --info=progress2 --remove-source-files ' tmp_dir '/* ' out_dir]);
end
  