%% Linear Models predicting ITNA Activity/Selectivity from connectivity
% Controlling for sex, age, and brain volume

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
[T,ind_subs_included,ind_subs_included_PP19] = innerjoin(T,PP19,'Keys','subID');

%% Specify models to run
              % Output Label            % Data File        
data_spec_map = {'FC_DigL_ALL_Zmap'     'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.ALL.niml.dset'
                 'FC_DigL_DpDa_Zdiff'   'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset'
                 'FC_DigL_LpLa_Zdiff'   'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Lp-La.niml.dset'
                 'FC_DigR_ALL_Zmap'     'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.ALL.niml.dset'
                 'FC_DigR_DpDa_Zdiff'   'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da.niml.dset'
                 'SC_DigL'              'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
                 'SC_LetL'              'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
                 'SC_DigR'              'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c.niml.dset'};        
%                  'SC_WhlBrn_Dens_lh'    'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.log.niml.dset'
%                  'SC_WhlBrn_Leng_lh'    'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset'
%                  'SC_WhlBrn_Dens_rh'    'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.log.niml.dset'
%                  'SC_WhlBrn_Leng_rh'    'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset'};

%% Run models
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
     P = ds.data(:,ind_subs_included); % Include only those subjects with data (130161 excluded from PP19 analyses due to motion)
     Stats1 = zeros(size(P,1),6); % = Model predicting Digit present activity vs baseline
     Stats2 = zeros(size(P,1),6); % = Model predicting Digit present > absent activity (i.e. selectivity)
     
     % Set the hemisphere of the data
     if contains(file_name,'lh')
         hemi = 'L';
     elseif contains(file_name,'rh')
         hemi = 'R';
     end
     
     % TRUE DATA 
     % Display a progress bar
     updateWaitbar = waitbarParfor(size(P,1), ['Calculation in progress for ' strrep(mdl_name,'_', ' ')]);
     % Parallel calculation of linear model for each node
     parfor jj = 1:size(P,1)
         Ttmp = T;
         Ttmp.Node_Conn = P(jj,:)';
         
         % Run linear modelling and save coefficients, tstats, and pvalues
         % Predict digit present activity vs baseline
         if contains(mdl_name,'SC')
             mdl_spec = [hemi '_ITG_D_Pres ~ Node_Conn + AgeMonths_T + female_T + BrainSegVolNotVent'];
         else
             mdl_spec = [hemi '_ITG_D_Pres ~ Node_Conn + AgeMonths_T + female_T'];
         end
         [est,t,p,bf10,bf01,bfboth] = run_linmdl(Ttmp,mdl_spec);
         Stats1(jj,:) = [est,t,p,bf10,bf01,bfboth];

         % Predict digit present > absent activity (selectivity)
         if contains(mdl_name,'SC')
             mdl_spec = [hemi '_ITG_CON ~ Node_Conn + AgeMonths_T + female_T + BrainSegVolNotVent'];
         else
             mdl_spec = [hemi '_ITG_CON ~ Node_Conn + AgeMonths_T + female_T'];
         end
         [est,t,p,bf10,bf01,bfboth] = run_linmdl(Ttmp,mdl_spec);
         Stats2(jj,:) = [est,t,p,bf10,bf01,bfboth];

         % Progress bar update
         updateWaitbar(); %#ok<PFBNS>
     end
     
     % Save stats surface niml file
     out_dir = [top_dir '/Results/LinMdl_' mdl_name];
     mkdir(out_dir);
     S = ds_mean;
     S.labels = {'Estimate','Tstat','Pval','BF10','BF01','BFBoth'};
     S.data = Stats1;
     afni_niml_writesimple(S,[out_dir '/LinMdl_Dp_x_' mdl_name '.niml.dset']);
     S.data = Stats2;
     afni_niml_writesimple(S,[out_dir '/LinMdl_DpDa_x_' mdl_name '.niml.dset']);
     
     
     
    % NULL DATA (based on shuffling residuals after nuisance regression)
    % *** Freedman and Lane procedure described Winkler et al. 2014.

    % Run linear NUISANCE modelling to get residuals
    % Predict digit present activity vs baseline
    if contains(mdl_name,'SC')
        mdl_spec = [hemi '_ITG_D_Pres ~ AgeMonths_T + female_T + BrainSegVolNotVent'];
    else
        mdl_spec = [hemi '_ITG_D_Pres ~ AgeMonths_T + female_T'];
    end
    L = fitlm(T,mdl_spec);
    T.Residuals_D_Pres = L.Residuals.Raw;
    % Predict digit present > absent activity (selectivity)
    if contains(mdl_name,'SC')
        mdl_spec = [hemi '_ITG_CON ~ AgeMonths_T + female_T + BrainSegVolNotVent'];
    else
        mdl_spec = [hemi '_ITG_CON ~ AgeMonths_T + female_T'];
    end
    L = fitlm(T,mdl_spec);
    T.Residuals_CON = L.Residuals.Raw;
    
    
    out_dir = [top_dir '/Results/LinMdl_' mdl_name];
    Stats1null = zeros(size(P,1),2); % = Model predicting Digit present activity vs baseline
    Stats2null = zeros(size(P,1),2); % = Model predicting Digit present > absent activity (i.e. selectivity)
    null_iters = 100;
    for nn = 1:null_iters
        % Randomize activity data
        Tnull = T;
        c = randperm(height(Tnull));
        Tnull.Residuals_D_Pres = Tnull.Residuals_D_Pres(c);
        Tnull.Residuals_CON = Tnull.Residuals_CON(c);
        
        %updateWaitbar = waitbarParfor(size(P,1), ['Calculation in progress for ' strrep(mdl_name,'_', ' ') ' - NULL iteration #' num2str(nn)]);
        parfor jj = 1:size(P,1)
            Ttmp = Tnull;
            Ttmp.Node_Conn = P(jj,:)';
            [b,SE] = linear_mdl_fast(Ttmp.Residuals_D_Pres,Ttmp.Node_Conn);
            Stats1null(jj,:) = [b(2), b(2)/SE(2)];
            [b,SE] = linear_mdl_fast(Ttmp.Residuals_CON,Ttmp.Node_Conn);
            Stats2null(jj,:) = [b(2), b(2)/SE(2)];
            % Progress bar update
            %updateWaitbar(); %#ok<PFBNS>
        end
        % Save data
        S = ds_mean;
        S.labels = {'Estimate','Tstat'};
        S.data = Stats1null;
        afni_niml_writesimple(S,[out_dir '/LinMdl_Dp_x_' mdl_name '_NULL_iter' num2str(nn, '%03.f') '.niml.dset']);
        S.data = Stats2null;
        afni_niml_writesimple(S,[out_dir '/LinMdl_DpDa_x_' mdl_name '_NULL_iter' num2str(nn, '%03.f') '.niml.dset']);
    end
end

%% Global linear modeling functions
function [est,t,p,bf10,bf01,bfboth] = run_linmdl(Ttmp,mdl_spec)
    L = fitlm(Ttmp,mdl_spec);
    ind = strcmp(L.CoefficientNames,'Node_Conn');
    est = L.Coefficients.Estimate(ind);
    t = L.Coefficients.tStat(ind);
    p = L.Coefficients.pValue(ind);
    % Get Bayes factors
    bf10 = bf.bfFromT(t,L.DFE);
    bf01 = 1/bf10;
    if bf01 >= 1
        bfboth = -bf01+1;
    elseif bf10 > 1
        bfboth = bf10-1;
    end   
end


function [b,SE] = linear_mdl_fast(Y,X)
% This function uses the fastest matrix operations to fit the model and get the
% coefficient standard errors (in order to compute tstats) 

%%Calculation of Standard Error With Intercept
n = length(X);                                      %  Number of observations
XBar=mean(X);                                       %  Calculates mean of X
YBar=mean(Y);                                       %  Calculates mean of Y
Sxx=sum((X-XBar).^2);  
Sxy=sum((X-XBar).*(Y-YBar));
Slope = Sxy/Sxx;                                    %  Calculates Slope
Intercept= YBar-Slope*XBar;                         %  Calculates Intercept
yfit=Intercept + X*Slope;                           %  Fitted response values based on the slope
r = Y - yfit;                                       %  r is the residuals, which is the observed minus fitted values
SSE = sum(r.^2);                                    %  SSE is the sum of squared errors
MSE=SSE/(n-2);                                      %  Mean Squared Error 
SE=[                                                %  Standard Error of the regression coefficients
    sqrt(MSE*sum(X.^2)/(n*Sxx));                    %  Standard Error of the intercept coefficient
    sqrt(MSE/Sxx)]  ;                                %  Standard Error of the slope coefficient
b = [Intercept,Slope];
%b =[ones(length(X),1),X]\Y;
end

