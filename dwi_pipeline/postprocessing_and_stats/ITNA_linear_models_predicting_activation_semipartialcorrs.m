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
PP19_Let = readtable('/Volumes/NBL_Projects/Price_NFA_RSA/PollackPrice2019/2_Data_and_Analyses/RSA/UnivariateAsymmetry/L_R_ITG_Let_betas_WCJ.csv');
PP19 = innerjoin(PP19,PP19_Let,'Keys','subject');
PP19.subID = string(PP19.subject);
[T,ind_subs_included,ind_subs_included_PP19] = innerjoin(T,PP19,'Keys','subID');

% Create double contrast variables
T.L_ITG_DBL_CON_PP19 = T.L_ITG_CON_PP19 - T.L_ITG_CON_PP19_Let;
T.R_ITG_DBL_CON_PP19 = T.R_ITG_CON_PP19 - T.R_ITG_CON_PP19_Let;

%% Specify models to run
              % Output Label            % Data File        
data_spec_map = {%'SC_DigL'              'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
                 %'SC_DigR'              'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c.niml.dset'
                 %'SC_LetL'              'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
                 'SC_DigL-DigR'         'AllSubs_tracks_ss3t_50M_Dp-Da_MINUS_Dp-Da_math.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'};    
             
%                  'FC_DigL_ALL_Zmap'     'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.ALL.niml.dset'
%                  'FC_DigL_DpDa_Zdiff'   'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset'
%                  'FC_DigL_LpLa_Zdiff'   'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Lp-La.niml.dset'
%                  'FC_DigR_ALL_Zmap'     'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.ALL.niml.dset'
%                  'FC_DigR_DpDa_Zdiff'   'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da.niml.dset'
%                  'SC_WhlBrn_Dens_lh'    'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.log.niml.dset'
%                  'SC_WhlBrn_Leng_lh'    'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset'
%                  'SC_WhlBrn_Dens_rh'    'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.log.niml.dset'
%                  'SC_WhlBrn_Leng_rh'    'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset'};

dependent_vars = {'calc_skills_ss_resid_PP19'
                  'calc_skills_ss_PP19'
                  'calcss_resid_PP19'
                  'lwidss_resid_PP19'
                  'fluencyss_resid_PP19'
                  'ITG_CON_PP19'     % Note the ITG variables will be appended with L_ or R_ depending on the input data
                  'ITG_DBL_CON_PP19'};
              
dependent_vars = {'LI_D_Pres'
                  'LI_CON_PP19'
                  'LI_D_Abs'};

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
     % Include connectivity values for only those subjects with data (130161 excluded from PP19 analyses due to motion)
     conn_data = ds.data(:,ind_subs_included); 
     
     % Loop through variables we are predicting (i.e., depedent variables)
     for pp = 1:numel(dependent_vars)
         % Populate variables
         dv = dependent_vars{pp};
         Stats = zeros(size(conn_data,1),6);
         
         % Set the hemisphere of the data
         if contains(file_name,'lh') && contains(dv,'ITG')
             dv = ['L_' dv];
         elseif contains(file_name,'rh') && contains(dv,'ITG')
             dv = ['R_' dv];
         end

         % TRUE DATA 
         % Display a progress bar
         updateWaitbar = waitbarParfor(size(conn_data,1), ['Calculation in progress for ' strrep(mdl_name,'_', ' ')]);
         % Parallel calculation of linear model for each node
         parfor jj = 1:size(conn_data,1)
             Ttmp = T;
             Ttmp.Node_Conn = conn_data(jj,:)';

%              % Run linear modelling and save coefficients, tstats, and pvalues
%              % Predict digit present activity vs baseline
%              if contains(mdl_name,'SC')
%                  mdl_spec = [dv ' ~ AgeMonths + female + BrainSegVolNotVent'];
%                  %mdl_spec = [dv ' ~ Node_Conn + AgeMonths + female + BrainSegVolNotVent'];
%              else
%                  mdl_spec = [dv ' ~ AgeMonths + female'];
%                  %mdl_spec = [dv ' ~ Node_Conn + AgeMonths + female'];
%              end
%              [est,t,p,bf10,bf01,bfboth] = run_linmdl(Ttmp,mdl_spec);
             
             % Run the semi-partial correlation
             [~,~,resid] = regress(Ttmp.Node_Conn, table2array(Ttmp(:,{'AgeMonths' 'female' 'BrainSegVolNotVent'})));  
             [bf10,r,p] = bf.corr(Ttmp.(dv),resid);
             n = height(Ttmp);
             t = r / ((1-r^2)/(n-2))^0.5
             bf01 = 1/bf10;
             if bf01 >= 1
                 bfboth = -bf01+1;
             elseif bf10 > 1
                 bfboth = bf10-1;
             end
             Stats(jj,:) = [r,t,p,bf10,bf01,bfboth];
             
             
%              if contains(mdl_name,'SC')
%                  X = [Ttmp.Node_Conn, Ttmp.AgeMonths , Ttmp.female, Ttmp.BrainSegVolNotVent];
%                  Y = Ttmp.(pv)
%              else
%                  X = [Ttmp.Node_Conn, Ttmp.AgeMonths, Ttmp.female];
%                  Y = Ttmp.(pv)
%              end
%              [b,t,p] = linear_mdl_fast(Y,X);
             
             % Progress bar update
             updateWaitbar(); %#ok<PFBNS>
         end

         % Save stats surface niml file
         out_dir = [top_dir '/Results/SemiPartCorr_' mdl_name];
         mkdir(out_dir);
         S = ds_mean;
         S.labels = {'Rho','Tstat','Pval','BF10','BF01','BFBoth'};
         S.data = Stats;
         afni_niml_writesimple(S,[out_dir '/SemiPartCorr_' dv '_x_' mdl_name '.niml.dset']);
         
        % NULL DATA (based on shuffling residuals after nuisance regression)
        % *** Freedman and Lane procedure described Winkler et al. 2014.

        % Run linear NUISANCE modelling to get residuals
        % Predict digit present activity vs baseline
%         if contains(mdl_name,'SC')
%             mdl_spec = [dv ' ~ AgeMonths + female + BrainSegVolNotVent'];
%         else
%             mdl_spec = [dv ' ~ AgeMonths + female'];
%         end
%         L = fitlm(T,mdl_spec);
%         T.Residuals_Tmp = L.Residuals.Raw;
    
        out_dir = [top_dir '/Results/SemiPartCorr_' mdl_name];
        StatsNull = zeros(size(conn_data,1),2); % = Model predicting Digit present activity vs baseline
        null_iters = 100;
        for nn = 1:null_iters
            % Randomize activity data
            Tnull = T;
            c = randperm(height(Tnull));
            %Tnull.Residuals_Tmp = Tnull.Residuals_Tmp(c);

            %updateWaitbar = waitbarParfor(size(P,1), ['Calculation in progress for ' strrep(mdl_name,'_', ' ') ' - NULL iteration #' num2str(nn)]);
            parfor jj = 1:size(conn_data,1)
                Ttmp = Tnull;
                Ttmp.Node_Conn = conn_data(jj,:)';
                
                % Run the semipartial correlation
                [~,~,resid] = regress(Ttmp.Node_Conn, table2array(Ttmp(:,{'AgeMonths' 'female' 'BrainSegVolNotVent'}))); 
                [r,p] = corr(Ttmp.(dv)(c),resid); % Randomize the order of the depedenent variable
                n = height(Ttmp);
                t = r / ((1-r^2)/(n-2))^0.5;
     
                %[b,SE] = linear_mdl_fast(Ttmp.Residuals_Tmp,Ttmp.Node_Conn);
                StatsNull(jj,:) = [r, t];
                % Progress bar update
                %updateWaitbar(); %#ok<PFBNS>
            end
            % Save data
            S = ds_mean;
            S.labels = {'Estimate','Tstat'};
            S.data = StatsNull;
            afni_niml_writesimple(S,[out_dir '/SemiPartCorr_' dv '_x_' mdl_name '_NULL_iter' num2str(nn, '%03.f') '.niml.dset']);
        end
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

%%Calculation of Coeffiecients and Standard Error With Intercept
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


% tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
% tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;              % 1-tailed t-distribution
% t = 0;
% p = 0;
end

