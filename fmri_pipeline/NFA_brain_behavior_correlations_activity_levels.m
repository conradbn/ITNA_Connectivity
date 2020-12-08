%% Brain-Behavior Correlations for ITNA Connectivity Project
% This script performs ALL the brain/behavior correlations for the paper,
% primarily involving the correlation of math ability with wholebrain
% activation to the digit condition, structural and functional connectivity
% of the digit ROIs, as well as with fiber density, fiber length, and
% activation within the ROIs. Nuisance variables are controlled for in some
% instances, including age, sex, and brain volume (when looking at
% structural metrics). Also, we are only considering connectivity of the
% ROIs within the ipsilateral hemispheres at this time. The following
% variables are considered (note that some are vectors "1D" while others
% are surface maps "2D"):
%
% 1)	Calculation Skills SS Residual ? 1D
% 2)	Digit Present Coefficient Digit L ? 1D
% 3)	Digit Present Coefficient Digit R ? 1D
% 4)	Digit Present ? Digit Absent Coefficient Digit L ? 1D
% 5)	Digit Present ? Digit Absent Coefficient Digit R ? 1D
% 6)	Mean streamline density Digit L ? 1D
% 7)	Mean streamline density Digit R ? 1D
% 8)	Mean streamline length Digit L ? 1D
% 9)	Mean streamline length Digit R ? 1D
% 10)	Brain volume ? 1D
% 11)	Sex ? 1D
% 12)	Age ? 1D
% 13)	Activation map Digit Present L ? 2D
% 14)	Activation map Digit Present R ? 2D
% 15)	Contrast map Digit Present ? Digit Absent L ? 2D
% 16)	Contrast map Digit Present ? Digit Absent R ? 2D
% 17)	Structural connectivity Digit L ? 2D
% 18)	Structural connectivity Digit R ? 2D
% 19)	Functional connectivity Zmap Digit L ? Digit Present ? 2D
% 20)	Functional connectivity Zmap Digit R ? Digit Present ? 2D
% 21)	Functional connectivity Zdiff Digit L ? Digit Present ? Digit Absent ? 2D
% 22)	Functional connectivity Zdiff Digit R ? Digit Present ? Digit Absent ? 2D


%% Setup
purge
top_dir = '/Users/benconrad/Desktop/BrainBehavCorrelations'; %'/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations';
cd(top_dir);
% task_conds = {'Dp' 'Da' 'Dp-Da'}; %'Da' 'Lp' 'La' 
% hemi = {'lh' 'rh'};

%% Create table with all 1D variables
% Behavioral & nuisance variables
T = readtable('Covariates/Export_Trimmed_wResidualizedScores_200619dy_check_DWI_29subj.csv');
T.Properties.VariableNames{'subject_nr'} = 'subID';

% Get the brain volume variable
Tvol = readtable('Covariates/aseg_stats_subid.csv');
Tvol = Tvol(:,{'subID','BrainSegVolNotVent'});
T = innerjoin(T,Tvol,'Keys','subID'); 

% Get the mean activation coefficients (and tstats)
task_conds = {'Dp' 'Da' 'Dp-Da'};%'Da' 'Lp' 'La' 
hemi = {'lh' 'rh'};


%% Load the ROI data
cd([top_dir '/AllSubs_dsets']);
roi_labels = {'Dig_lh_141' 'AllSubs_std.141.lh.PollackandPrice19_Dp-Da.MNI152.votc.inflated.14mm_diam.niml.dset'
              'Let_lh_141' 'AllSubs_std.141.lh.PollackandPrice19_Dp-Da.MNI152.votc.inflated.14mm_diam.niml.dset'
              'Dig_rh_141' 'AllSubs_std.141.rh.PP19_Dp-Da_math.MNI152.votc.inflated.14mm_diam.niml.dset'
              'Dig_lh_60' 'AllSubs_std.60.lh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.niml.dset'
              'Let_lh_60' 'AllSubs_std.60.lh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam.niml.dset'
              'Dig_rh_60' 'AllSubs_std.60.rh.PP19_Dp-Da_math.MNI152.votc.inflated.14mm_diam.niml.dset'};
    
for ii = 1:size(roi_labels,1)
    r = afni_niml_readsimple(roi_labels{ii,2});
    rois.(roi_labels{ii,1}) = r.data;
end

%% Extract mean data from ROIs
% Specify data to extract
              % Output Label         % ROI Label     % Data File        
data_spec_roi = {'Coeff_Dp_Da_Dig_rh'   'Dig_rh_60'     'AllSubs_stats_rh_REML_Coeff_Dp-Da.niml.dset'
                 'Coeff_Dp_Da_Dig_lh'   'Dig_lh_60'     'AllSubs_stats_lh_REML_Coeff_Dp-Da.niml.dset'
                 'Tstat_Dp_Da_Dig_rh'   'Dig_rh_60'     'AllSubs_stats_rh_REML_Tstat_Dp-Da.niml.dset'
                 'Tstat_Dp_Da_Dig_lh'   'Dig_lh_60'     'AllSubs_stats_lh_REML_Tstat_Dp-Da.niml.dset'
                 'Coeff_Lp_La_Dig_rh'   'Dig_rh_60'     'AllSubs_stats_rh_REML_Coeff_Lp-La.niml.dset'
                 'Coeff_Lp_La_Dig_lh'   'Dig_lh_60'     'AllSubs_stats_lh_REML_Coeff_Lp-La.niml.dset'
                 'Tstat_Lp_La_Dig_rh'   'Dig_rh_60'     'AllSubs_stats_rh_REML_Tstat_Lp-La.niml.dset'
                 'Tstat_Lp_La_Dig_lh'   'Dig_lh_60'     'AllSubs_stats_lh_REML_Tstat_Lp-La.niml.dset'
                 'Coeff_Dp_Dig_rh'      'Dig_rh_60'     'AllSubs_stats_rh_REML_Coeff_Dp.niml.dset'
                 'Coeff_Dp_Dig_lh'      'Dig_lh_60'     'AllSubs_stats_lh_REML_Coeff_Dp.niml.dset'
                 'Tstat_Dp_Dig_rh'      'Dig_rh_60'     'AllSubs_stats_rh_REML_Tstat_Dp.niml.dset'
                 'Tstat_Dp_Dig_lh'      'Dig_lh_60'     'AllSubs_stats_lh_REML_Tstat_Dp.niml.dset'
                 'Coeff_Da_Dig_rh'      'Dig_rh_60'     'AllSubs_stats_rh_REML_Coeff_Da.niml.dset'
                 'Coeff_Da_Dig_lh'      'Dig_lh_60'     'AllSubs_stats_lh_REML_Coeff_Da.niml.dset'
                 'Tstat_Da_Dig_rh'      'Dig_rh_60'     'AllSubs_stats_rh_REML_Tstat_Da.niml.dset'
                 'Tstat_Da_Dig_lh'      'Dig_lh_60'     'AllSubs_stats_lh_REML_Tstat_Da.niml.dset'
                 'Coeff_Lp_Dig_rh'      'Dig_rh_60'     'AllSubs_stats_rh_REML_Coeff_Lp.niml.dset'
                 'Coeff_Lp_Dig_lh'      'Dig_lh_60'     'AllSubs_stats_lh_REML_Coeff_Lp.niml.dset'
                 'Tstat_Lp_Dig_rh'      'Dig_rh_60'     'AllSubs_stats_rh_REML_Tstat_Lp.niml.dset'
                 'Tstat_Lp_Dig_lh'      'Dig_lh_60'     'AllSubs_stats_lh_REML_Tstat_Lp.niml.dset'
                 'Coeff_La_Dig_rh'      'Dig_rh_60'     'AllSubs_stats_rh_REML_Coeff_La.niml.dset'
                 'Coeff_La_Dig_lh'      'Dig_lh_60'     'AllSubs_stats_lh_REML_Coeff_La.niml.dset'
                 'Tstat_La_Dig_rh'      'Dig_rh_60'     'AllSubs_stats_rh_REML_Tstat_La.niml.dset'
                 'Tstat_La_Dig_lh'      'Dig_lh_60'     'AllSubs_stats_lh_REML_Tstat_La.niml.dset'
                 'Dens_Dig_lh'          'Dig_lh_141'    'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.log.niml.dset'
                 'Leng_Dig_lh'          'Dig_lh_141'    'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset'
                 'Dens_Dig_rh'          'Dig_rh_141'    'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.log.niml.dset'
                 'Leng_Dig_rh'          'Dig_rh_141'    'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset'};

 % Loop through each dataset specified (row)
 for ii = 1:size(data_spec_roi,1)
     % Get labels and data file
     l = data_spec_roi{ii,1};
     r = logical(rois.(data_spec_roi{ii,2}));
     ds = afni_niml_readsimple(data_spec_roi{ii,3});
     d = ds.data;
     m = [];
     for jj = 1:size(d,2)
         m(jj,1) = mean(d(r(:,jj),jj));
     end
     % Add to table
     T.(l) = m;
 end

%% Load the 2D surface data (i.e. activity and connectivity maps)
% Specify data to extract
              % Output Label                    % Data File        
data_spec_map = {'Coeff_Dp_Da_Dig_rh_2D'       'AllSubs_stats_rh_REML_Coeff_Dp-Da.niml.dset'
                 'Coeff_Dp_Da_Dig_lh_2D'       'AllSubs_stats_lh_REML_Coeff_Dp-Da.niml.dset'
                 'Tstat_Dp_Da_Dig_rh_2D'       'AllSubs_stats_rh_REML_Tstat_Dp-Da.niml.dset'
                 'Tstat_Dp_Da_Dig_lh_2D'       'AllSubs_stats_lh_REML_Tstat_Dp-Da.niml.dset'
                 'Dens_lh_2D'                  'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.log.niml.dset'
                 'Leng_lh_2D'                  'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset'
                 'Dens_rh_2D'                  'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.log.niml.dset'
                 'Leng_rh_2D'                  'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset'
                 'FConn_Dig_Zdiff_Dp_Da_lh_2D' 'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset'
                 'FConn_Dig_Zdiff_Dp_Da_rh_2D' 'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da.niml.dset'
                 'FConn_Dig_Zmap_Dp_lh_2D'     'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset'
                 'FConn_Dig_Zmap_Dp_rh_2D'     'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp.niml.dset'      
                 'SConn_Dig_lh_2D'             'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
                 'SConn_Dig_rh_2D'             'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c.niml.dset'};
clear Maps             
% Loop through each dataset specified (row)
for ii = 1:size(data_spec_map,1)
     % Get labels and data file
     l = data_spec_map{ii,1};
     ds = afni_niml_readsimple(data_spec_map{ii,2});
     % Add to structure variable
     Maps.(l) = ds.data;
end

%% Run Vector-Map Brain-Behavior Correlations of Interest
% % Specify data to extract
%                   % Vector         % Map                         % NuisanceVars
% data_spec_corr = {'calc_skills_ss' 'SConn_Dig_lh_2D'             'AgeMonths,female,BrainSegVolNotVent' 
%                   'calc_skills_ss' 'SConn_Dig_rh_2D'             'AgeMonths,female,BrainSegVolNotVent'
%                   'calc_skills_ss' 'FConn_Dig_Zmap_Dp_lh_2D'     'AgeMonths,female' 
%                   'calc_skills_ss' 'FConn_Dig_Zdiff_Dp_Da_lh_2D' 'AgeMonths,female'};
%               
%               
% % % Loop through each dataset specified (row)
% % for ii = 1:size(data_spec_corr,1)
% %      % Get labels and data file
% %      l = [data_spec_corr{ii,1} '_X_' data_spec_corr{ii,2}];
% %      v = T.(data_spec_corr{ii,1});
% %      d = Maps.(data_spec_corr{ii,2});
% %      cv = data_spec_corr{ii,1};
% %      for
% %      % Calculate correlation map
% %      C = corr(v,d');
% % end

%% Seed-seed Structural connectivity matrices and contrast of LH vs RH
Mats.SConn_Dig_lh_1D = readmatrix('AllSubs_tracks_ss3t_50M_Dp-Da.lh.fingerprint.norm.csv');
Mats.SConn_Dig_rh_1D = readmatrix('AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.fingerprint.norm.csv');
lut = readmatrix('aparc.a2009s+aseg_REN_all.niml.lt','FileType','text','Range','6:202','OutputType','char');

% Get unique region names
ctx = lut(contains(lut(:,2),'ctx_'),2);
ctx_rois = cellfun(@(S) S(8:end), ctx, 'Uniform', 0);
ctx_rois = unique(ctx_rois);
ctx_rois(contains(ctx_rois,'Unknown')) = [];
Mats.ROIs = ctx_rois;

% Get SUMA MNI FreeSurfer Parcellation information
FS = afni_niml_read('../FreeSurfer_ROIs/std.141.lh.aparc.a2009s.annot.niml.dset');
FS_rois = FS.nodes{1}.data;
FS_nodes = FS.nodes{2}.data;
ctx_rois_ids = readmatrix('../FreeSurfer_ROIs/FS_idcodes.csv','OutputType','char');
Tval_map = afni_niml_readsimple('../GroupMean_dsets/GroupMean_std.141.lh.PollackandPrice19_Lp-La.MNI152.votc.inflated.14mm_diam.niml.dset');
Tval_map.data(:) = 0; 

% Loop through each ROI, run paired t-tests
for ii = 1:numel(ctx_rois)
    c = ctx_rois{ii};
    rh = ['ctx_rh_' c];
    lh = ['ctx_lh_' c];
    rh_ind = str2double(lut{strcmp(lut(:,2),rh),1});
    lh_ind = str2double(lut{strcmp(lut(:,2),lh),1});
    % Paired t-test on connectivity values
    [h,p,ci,stats] = ttest(Mats.SConn_Dig_lh_1D(lh_ind,:),Mats.SConn_Dig_rh_1D(rh_ind,:));
    % Set the ROI ID code for plotting
    id = ctx_rois_ids(contains(ctx_rois_ids(:,2),['lh_' c]),1);
    % Associate T-values with surface nodes
    Tval_map.data(FS_rois == str2double(id)) = stats.tstat;
    % Save the info in structure
    Mats.Tstats(ii,1) = stats.tstat;
    Mats.Pvals(ii,1) = p;
    Mats.ROI_IDs(ii,1) = id;
end


% Write Tvalue map
afni_niml_writesimple(Tval_map,'../Results/SConn_lh_vs_rh_Tvals.niml.dset')

% BH FDR correction
Mats.Pvals_FDR = mafdr(Mats.Pvals,'BHFDR',true);

% Write map only significnat Tvalues
Tval_map_sig_05 = Tval_map;
Tval_map_sig_05.data(:) = 0;
Tval_map_sig_01 = Tval_map;
Tval_map_sig_01.data(:) = 0;

for ii = 1:numel(ctx_rois)
    p = Mats.Pvals_FDR(ii);
    t = Mats.Tstats(ii);
    id = Mats.ROI_IDs(ii);
    if  p < 0.05
        Tval_map_sig_05.data(FS_rois == str2double(id)) = t;
    end
    if p < 0.01
        Tval_map_sig_01.data(FS_rois == str2double(id)) = t;
    end
end

afni_niml_writesimple(Tval_map_sig_05,'../Results/SConn_lh_vs_rh_Tvals_FDR_p05.niml.dset')
afni_niml_writesimple(Tval_map_sig_01,'../Results/SConn_lh_vs_rh_Tvals_FDR_p01.niml.dset')


%% Seed-seed Functional connectivity matrices and contrast of LH vs RH
FS_lh = afni_niml_read('../FreeSurfer_ROIs/std.60.lh.aparc.a2009s.annot.niml.dset');
FS_lh_rois = FS_lh.nodes{1}.data;
FS_rh = afni_niml_read('../FreeSurfer_ROIs/std.60.lh.aparc.a2009s.annot.niml.dset');
FS_rh_rois = FS_rh.nodes{1}.data;

lh_conn = Maps.FConn_Dig_Zmap_Dp_lh_2D; %Maps.FConn_Dig_Zdiff_Dp_Da_lh_2D
rh_conn = Maps.FConn_Dig_Zmap_Dp_rh_2D; %Maps.FConn_Dig_Zdiff_Dp_Da_rh_2D

% USE THE SAME 141 DENSITY MAP FROM ABOVE
%Tval_map = afni_niml_readsimple('../GroupMean_dsets/GroupMean_stats_lh_REML_Tstat_Dp.niml.dset');
Tval_map.data(:) = 0; 

for ii = 1:numel(Mats.ROI_IDs)
    id = str2double(Mats.ROI_IDs(ii));
    
    ind_lh = FS_lh_rois == id;
    ind_rh = FS_rh_rois == id;
    
    lh_conn_mean = mean(lh_conn(ind_lh,:),1)';
    rh_conn_mean = mean(rh_conn(ind_rh,:),1)';
    
    % Ttests
    [h,p,ci,stats] = ttest(lh_conn_mean,rh_conn_mean);
    
    Tval_map.data(FS_rois == id) = stats.tstat;
    
    % Save to structure
    Mats.FConn_Dig_lh_Dp_1D(ii,:) = lh_conn_mean;
    Mats.FConn_Dig_rh_Dp_1D(ii,:) = rh_conn_mean;
    Mats.FConn_Tstats(ii,1) = stats.tstat;
    Mats.FConn_Pvals(ii,1) = p;
    
end


% Write Tvalue map
afni_niml_writesimple(Tval_map,'../Results/FConn_lh_vs_rh_Dp_Tvals.niml.dset')

% BH FDR correction
Mats.FConn_Pvals_FDR = mafdr(Mats.FConn_Pvals,'BHFDR',true);

% Write map only significnat Tvalues
Tval_map_sig_05 = Tval_map;
Tval_map_sig_05.data(:) = 0;
Tval_map_sig_01 = Tval_map;
Tval_map_sig_01.data(:) = 0;

for ii = 1:numel(ctx_rois)
    p = Mats.FConn_Pvals_FDR(ii);
    t = Mats.FConn_Tstats(ii);
    id = Mats.ROI_IDs(ii);
    if  p < 0.05
        Tval_map_sig_05.data(FS_rois == str2double(id)) = t;
    end
    if p < 0.01
        Tval_map_sig_01.data(FS_rois == str2double(id)) = t;
    end
end

afni_niml_writesimple(Tval_map_sig_05,'../Results/FConn_lh_vs_rh_Dp_Tvals_FDR_p05.niml.dset')
afni_niml_writesimple(Tval_map_sig_01,'../Results/FConn_lh_vs_rh_Dp_Tvals_FDR_p01.niml.dset')

%% Linear modeling predicting Digit response in each ITNA
% Separate models for Structural and Functional connectivity
% Covariates - Age, Sex, Brain Volume, Lp-La
                % Table Variable           % Example file from which to adopt structure
Predictor_Map = {'FConn_Dig_Zmap_Dp_lh_2D'     '../GroupMean_dsets/GroupMean_stats_lh_REML_Tstat_Dp.niml.dset'
                 'FConn_Dig_Zmap_Dp_rh_2D'     '../GroupMean_dsets/GroupMean_stats_rh_REML_Tstat_Dp.niml.dset'
                 'FConn_Dig_Zdiff_Dp_Da_lh_2D' '../GroupMean_dsets/GroupMean_stats_lh_REML_Tstat_Dp.niml.dset'
                 'FConn_Dig_Zdiff_Dp_Da_rh_2D' '../GroupMean_dsets/GroupMean_stats_rh_REML_Tstat_Dp.niml.dset'
                 'SConn_Dig_lh_2D'             '../GroupMean_dsets/GroupMean_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
                 'SConn_Dig_rh_2D'             '../GroupMean_dsets/GroupMean_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'};
             
for ii = 1:size(Predictor_Map,1)
    disp(num2str(ii));
    P = Maps.(Predictor_Map{ii,1});
    StatsC = [];
    StatsT = [];
    StatsP = [];
    parfor jj = 1:size(P,1)
        Ttmp = T;
        Ttmp.Node_Conn = P(jj,:)';
        % Run linear modelling and save coefficients, tstats, and pvalues
        if contains(Predictor_Map{ii,1},'lh')
            L = fitlm(Ttmp,'Coeff_Dp_Da_Dig_lh ~ Node_Conn + AgeMonths + female + BrainSegVolNotVent');
            ind = strcmp(L.CoefficientNames,'Node_Conn');
            StatsC(jj) = L.Coefficients.Estimate(ind);
            StatsT(jj) = L.Coefficients.tStat(ind);
            StatsP(jj) = L.Coefficients.pValue(ind);
        elseif contains(Predictor_Map{ii,1},'rh')
            L = fitlm(Ttmp,'Coeff_Dp_Da_Dig_rh ~ Node_Conn + AgeMonths + female + BrainSegVolNotVent');
            ind = strcmp(L.CoefficientNames,'Node_Conn');
            StatsC(jj) = L.Coefficients.Estimate(ind);
            StatsT(jj) = L.Coefficients.tStat(ind);
            StatsP(jj) = L.Coefficients.pValue(ind);
        end
    end
    StatsOut = [];
    StatsOut(:,1) = StatsC';
    StatsOut(:,2) = StatsT';
    StatsOut(:,3) = StatsP';
    % Save to surface niml file
    S = afni_niml_readsimple(Predictor_Map{ii,2});
    S.data = StatsOut;
    afni_niml_writesimple(S,['../Results/LinMdl_Dp-Da_Coeff_NOLETTERCONTROL_x_' Predictor_Map{ii,1} '.niml.dset']);
end







