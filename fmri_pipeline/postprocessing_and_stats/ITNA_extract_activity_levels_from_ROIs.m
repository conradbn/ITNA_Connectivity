%% Extract the beta weights (activity levels) for all subjects
purge
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData';
cd(start_dir);
sub_dirs = dir('*_proc');
task_conds = {'Dp' 'Da' 'Lp' 'La' 'Dp-Da' 'Lp-La'};
seed = {'Dp-Da' 'Lp-La'};
hemi = {'lh' 'rh'};

% Populate empty table
T = table();

% For each subject
for ii = 1:numel(sub_dirs)
    disp(['Working on subject #' num2str(ii) ' of ' num2str(numel(sub_dirs))]);
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
    

    %% Extract stats data (beta coefficients)
    % For each hemisphere
    for hh = 1:numel(hemi)
        h = hemi{hh};
        % Read in the stats dataset
        stats = [dir_results 'stats.' subj '.' h '_REML.niml.dset'];
        stats = afni_niml_readsimple(stats);

        % For each seed ROI
        for ss = 1:numel(seed)
            s = seed{ss};
            % Get the seed (ROI) mask node indices based on hemisphere/condition
            sm = [dir_fs 'SUMA/std.60.' h '.PP19_' s '.MNI152.votc.inflated.14mm_diam.1.1D'];
            sm = readmatrix(sm,'NumHeaderLines',2,'FileType','text');
            seed_inds = sm(:,1)+1;
            
            % For each task condition/contrast
            for tt = 1:numel(task_conds)
                t = task_conds{tt};
                brikT = strcmp(stats.labels,[t '#0_Tstat']);
                brikB = strcmp(stats.labels,[t '#0_Coef']);
                varnameT = [s '_seed_' t '_data_' h '_tstat'];
                varnameB = [s '_seed_' t '_data_' h '_coef'];
                T.(varnameT)(ii) = mean(stats.data(seed_inds,brikT));
                T.(varnameB)(ii) = mean(stats.data(seed_inds,brikB));
            end
        end
    end
end


%% Correlations with math
close all
T_activations = T;
covars = readtable(['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Covariates/Export_Trimmed_wResidualizedScores_200619dy_check_DWI_29subj.csv']);  
covars_dwi = load('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/TDI_endpoint_density_and_length.mat');
covars_dwi = covars_dwi.T;
covars_darren = readtable('/Volumes/NBL_Projects/Price_NFA_RSA/PollackPrice2019/2_Data_and_Analyses/RSA/UnivariateAsymmetry/L_R_ITG_Dig_betas_WCJ.csv');
rm = [1,2,3,30];
covars_darren(rm,:) = [];

rm = [10];
covars(rm,:) = [];
covars_dwi(rm,:) = [];
T_activations(rm,:) = [];
    
[R,P] = corr(covars.calc_skills_ss_resid,T_activations.("Dp-Da_seed_Dp-Da_data_lh_coef"));
[R,P] = corr(covars.calc_skills_ss_resid,T_activations.("Dp-Da_seed_Dp-Da_data_lh_tstat"));

[R,P] = corrcoef([covars.calc_skills_ss_resid,covars_dwi.Variables]);

notBoxPlot([T_activations.("Dp-Da_seed_Dp_data_lh_coef"),T_activations.("Dp-Da_seed_Da_data_lh_coef")...
            T_activations.("Dp-Da_seed_Dp_data_rh_coef"),T_activations.("Dp-Da_seed_Da_data_rh_coef")]);

        
        
[R,P] = corr(covars.calc_skills_ss_resid,covars_darren.calc_skills_ss_resid)

[R,P] = corr(T_activations.("Dp-Da_seed_Dp-Da_data_lh_coef"),covars_darren.L_ITG_CON)
[R,P] = corr(T_activations.("Dp-Da_seed_Dp-Da_data_rh_coef"),covars_darren.R_ITG_CON)
[R,P] = corr(covars.calc_skills_ss_resid,T_activations.("Dp-Da_seed_Dp-Da_data_rh_coef"))  
[R,P] = corr(covars.calc_skills_ss_resid,covars_darren.R_ITG_CON) 



[R,P] = corr(T_activations.("Dp-Da_seed_Dp-Da_data_rh_coef"),covars_darren.R_ITG_CON)

      
        






