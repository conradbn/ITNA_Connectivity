%% Extract the beta weights (activity levels) for all subjects
purge
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/GroupLevel/BrainBehavCorrelations';
cd(start_dir);
task_conds = {'Dp' 'Da' 'Dp-Da'};%'Da' 'Lp' 'La' 
hemi = {'lh' 'rh'};

covar = {'calc_skills_ss_resid','calc_skills_ss'};
covars = readtable(['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Covariates/Export_Trimmed_wResidualizedScores_200619dy_check_DWI_29subj.csv']);  


for hh = 1:numel(hemi)
    h = hemi{hh};
    unix(['3dTcat -overwrite -prefix all_sub_activation_stats.' h '.niml.dset /Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/*_proc/*.results/stats.1*.' h '_REML.niml.dset']);
    stats = afni_niml_readsimple(['all_sub_activation_stats.' h '.niml.dset']);
    for tt = 1:numel(task_conds)
        t = task_conds{tt};
        brikT = strcmp(stats.labels,[t '#0_Tstat']);
        brikB = strcmp(stats.labels,[t '#0_Coef']);
        statsT = stats.data(:,brikT);
        statsB = stats.data(:,brikB);
        for cc = 1:numel(covar)
            c = covar{cc};
            c_vec = covars.(c);
            % Run correlation - with Activation Coefficients
            corr_data = corr(statsB',c_vec);
            corr_data_Z = atanh(corr_data);
            % Write data
            out = stats;
            out.data = corr_data_Z;
            out_name = ['corr_FishZ.' c '.x.' t '.' h '_Coef.niml.dset'];
            afni_niml_writesimple(out,out_name);
%% Correlations with math
% close all
            
            % Run correlation - with Activation Tstats
            corr_data = corr(statsT',c_vec);
            corr_data_Z = atanh(corr_data);
            % Write data
            out = stats;
            out.data = corr_data_Z;
            out_name = ['corr_FishZ.' c '.x.' t '.' h '_Tstat.niml.dset'];
            afni_niml_writesimple(out,out_name);
        end
    end
end




% T_activations = T;
% covars = readtable(['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Covariates/Export_Trimmed_wResidualizedScores_200619dy_check_DWI_29subj.csv']);  
% covars_dwi = load('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/TDI_endpoint_density_and_length.mat');
% covars_dwi = covars_dwi.T;
% covars_darren = readtable('/Volumes/NBL_Projects/Price_NFA_RSA/PollackPrice2019/2_Data_and_Analyses/RSA/UnivariateAsymmetry/L_R_ITG_Dig_betas_WCJ.csv');
% rm = [1,2,3,30];
% covars_darren(rm,:) = [];
% 
% rm = [10];
% covars(rm,:) = [];
% covars_dwi(rm,:) = [];
% T_activations(rm,:) = [];
%     
% [R,P] = corr(covars.calc_skills_ss_resid,T_activations.("Dp-Da_seed_Dp-Da_data_lh_coef"));
% [R,P] = corr(covars.calc_skills_ss_resid,T_activations.("Dp-Da_seed_Dp-Da_data_lh_tstat"));
% 
% [R,P] = corrcoef([covars.calc_skills_ss_resid,covars_dwi.Variables]);
% 
% notBoxPlot([T_activations.("Dp-Da_seed_Dp_data_lh_coef"),T_activations.("Dp-Da_seed_Da_data_lh_coef")...
%             T_activations.("Dp-Da_seed_Dp_data_rh_coef"),T_activations.("Dp-Da_seed_Da_data_rh_coef")]);
% 
%         
%         
% [R,P] = corr(covars.calc_skills_ss_resid,covars_darren.calc_skills_ss_resid)
% 
% [R,P] = corr(T_activations.("Dp-Da_seed_Dp-Da_data_lh_coef"),covars_darren.L_ITG_CON)
% [R,P] = corr(T_activations.("Dp-Da_seed_Dp-Da_data_rh_coef"),covars_darren.R_ITG_CON)
% [R,P] = corr(covars.calc_skills_ss_resid,T_activations.("Dp-Da_seed_Dp-Da_data_rh_coef"))  
% [R,P] = corr(covars.calc_skills_ss_resid,covars_darren.R_ITG_CON) 
% 
% 
% 
% [R,P] = corr(T_activations.("Dp-Da_seed_Dp-Da_data_rh_coef"),covars_darren.R_ITG_CON)
% 
%       
%         
% 
% 
% 
% 
% 
% 
