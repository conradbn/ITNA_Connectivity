
%% Setup
purge
top_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper'; 
cd(top_dir);

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

%% Merge activation stats from Pollack & Price 2019 into table
% Read in activation stats from Darren's BV analysis
PP19 = readtable('/Volumes/NBL_Projects/Price_NFA_RSA/PollackPrice2019/2_Data_and_Analyses/RSA/UnivariateAsymmetry/L_R_ITG_Dig_betas_WCJ.csv');
PP19.subID = string(PP19.subject);
T = innerjoin(T,PP19,'Keys','subID');

%% Linear Models predicting ITNA Activity/Selectivity from connectivity
% Specify data to extract
              % Output Label            % Data File        
data_spec_map = {'SC_WhlBrn_Dens_lh'    'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.log.niml.dset'
                 'SC_WhlBrn_Leng_lh'    'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset'
                 'SC_WhlBrn_Dens_rh'    'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.log.niml.dset'
                 'SC_WhlBrn_Leng_rh'    'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset'
                                      
                 'SC_DigL'              'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
                 'SC_LetL'              'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
                 'SC_DigR'              'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c.niml.dset'
                 
                 'FC_DigL_ALL_Zmap'     'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.ALL.niml.dset'
                 'FC_DigL_DpDa_Zdiff'   'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset'
                 'FC_LetL_ALL_Zmap'     'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.ALL.niml.dset'
                 'FC_LetL_LpLa_Zdiff'   'AllSubs_Lp-La.lh.beta_series_corr.lh.Zdiff.Lp-La.niml.dset'
                 'FC_DigR_ALL_Zmap'     'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.ALL.niml.dset'
                 'FC_DigR_DpDa_Zdiff'   'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da.niml.dset'};
             
clear Maps             
% Loop through each dataset specified (row) and run linear model predicting
% activity levels
for ii = 1:size(data_spec_map,1)
     disp(num2str(ii));
     % Get labels and data file, as well as group mean for output data
     % template
     mdl_name = data_spec_map{ii,1};
     file_name = data_spec_map{ii,2};
     ds = afni_niml_readsimple([top_dir '/AllSubs_dsets/' file_name]);
     ds_mean = afni_niml_readsimple([top_dir '/GroupMean_dsets/' strrep(file_name,'AllSubs','GroupMean')]);
     % Populate variables
     P = ds.data;
     Stats1 = zeros(size(P,1),3);
     Stats2 = zeros(size(P,1),3);
     parfor jj = 1:size(P,1)
         Ttmp = T;
         Ttmp.Node_Conn = P(jj,:)';
         % Run linear modelling and save coefficients, tstats, and pvalues
         if contains(file_name,'lh')
             L = fitlm(Ttmp,'L_ITG_D_Pres ~ Node_Conn + AgeMonths_T + female_T + BrainSegVolNotVent');
             ind = strcmp(L.CoefficientNames,'Node_Conn');
             Stats1(jj,:) = [L.Coefficients.Estimate(ind), L.Coefficients.tStat(ind), L.Coefficients.pValue(ind)];
             L = fitlm(Ttmp,'L_ITG_CON ~ Node_Conn + AgeMonths_T + female_T + BrainSegVolNotVent');
             ind = strcmp(L.CoefficientNames,'Node_Conn');
             Stats2(jj,:) = [L.Coefficients.Estimate(ind), L.Coefficients.tStat(ind), L.Coefficients.pValue(ind)];
         elseif contains(file_name,'rh')
             L = fitlm(Ttmp,'R_ITG_D_Pres ~ Node_Conn + AgeMonths_T + female_T + BrainSegVolNotVent');
             ind = strcmp(L.CoefficientNames,'Node_Conn');
             Stats1(jj,:) = [L.Coefficients.Estimate(ind), L.Coefficients.tStat(ind), L.Coefficients.pValue(ind)];
             L = fitlm(Ttmp,'R_ITG_CON ~ Node_Conn + AgeMonths_T + female_T + BrainSegVolNotVent');
             ind = strcmp(L.CoefficientNames,'Node_Conn');
             Stats2(jj,:) = [L.Coefficients.Estimate(ind), L.Coefficients.tStat(ind), L.Coefficients.pValue(ind)];
         end
     end
     % Save to surface niml file
     out_dir = [top_dir '/Results/LinMdl_' mdl_name];
     mkdir(out_dir);
     S = ds_mean;
     S.data = StatsOut1;
     afni_niml_writesimple(S,[out_dir '/LinMdl_Dp_x_' mdl_name '.niml.dset']);
     S = ds_mean;
     S.data = StatsOut2;
     afni_niml_writesimple(S,[out_dir '/LinMdl_DpDa_x_' mdl_name '.niml.dset']);
end