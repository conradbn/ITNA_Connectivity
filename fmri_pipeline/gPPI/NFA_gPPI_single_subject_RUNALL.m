%% Run the gPPI analysis pipeline for all subjects
purge
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData';
cd(start_dir);
sub_dirs = dir('*_proc');

% For each subject
for ii = 1:numel(sub_dirs)
    cd(start_dir)
    sd = sub_dirs(ii).name;
    
    % Go to subject directory and load required variables
    cd(sd)
    var_file = dir('workspace_variables_post*.mat');
    load(var_file.name,'subj','dir_results','dir_nii','dir_fs');
    
    % Correct folder paths to direct to server
    dir_results = strrep(dir_results,'/Users/nbl_imac/Documents/afni_tmp/','/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/');
    dir_nii = strrep(dir_nii,'/Users/nbl_imac/Documents/afni_tmp/','/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/');
    dir_fs = strrep(dir_fs,'/Users/nbl_imac/Documents/afni_tmp/','/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/');
    
    % Run the gPPI pipeline
    NFA_gPPI_single_subject(subj,dir_nii,dir_results,dir_fs)
end