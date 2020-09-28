%% Setup
% Set subject
subj = '130843';
% Set number of runs
nruns = 1:5;
% Set number of time points per run in TR
n_tp = 480;
TR = 2;

% Up-sample the data because of stimulus duration of 3s
sub_TR = 1;
% Set seed label
sd = 'amygdala';
% Set condition list three conditions
condList = {'A' 'B' 'C'};

%% Create Gamma impulse response function
unix(['waver -dt ' sub_TR ' -GAM -peak 1 -inline 1@1 > GammaHR.1D']);

%% For each run, extract seed time series, run deconvolution, and create interaction regressor
for ii = 1:numel(nruns)
    r = nruns(ii); % foreach cc (`count -digits 1 1 $nruns`)
    % 1. For each run extract the average time series of the ROI (if afni_proc.py was used before, consider using the pb04.* file from the pipeline as input here)
    % 2. Remove the trend from the seed time series
    % The output SeedR.1D is a one-row text file. Convert the one-row time series to one column:
    % 2a. If your stimulus onset times were not synchronized with TR grids, pick up a sub_TR, e.g., 0.1 seconds, replace the above 1dtranspose step and upsample seed time series by xx (original TR divided by sub_TR) times:
    % 3. Run deconvolution on the seed time series
    % First generate the impulse response function, then run 3dTfitter
    
    
    


    unix(['3dmaskave -mask ROI+orig -quiet pb04.$subj.r0' r '.scale+tlrc > Seed' r '_' sd '.1D'...
          ' 3dDetrend -polort 2 -prefix SeedR' r '_' sd ' Seed' r '_' sd '.1D'...
          ' 1dtranspose SeedR' r '_' sd '.1D Seed_ts' r '_' sd 'D.1D'...
          ' rm -f SeedR' r '_' sd '.1D'...
          ' 1dUpsample 2 Seed_ts' r '_' sd 'D.1D > Seed_ts' r '_' sd '.1D'...
          ' 3dTfitter -RHS Seed_ts' r '_' sd '.1D -FALTUNG GammaHR.1D'...
          ' Seed_Neur' r '_' sd ' 012 -1']);
      for jj = 1:numel(condList)
          c = condList(jj);
          % First create a 1D (one column) file, conditionA.1D,  with 0's
          % (at those TR's where condition A does not occur), 1's (at those
          % TR's where condition A occurs). Unlike the traditional PPI, no
          % contrasting is performed here (thus no -1's).
          % The interaction is then created using waver
          
          unix(['head -' r ' stimuli/Allruns_stim_' c '_time.1D |tail -1 > tmp.1D'...
                'waver -dt ' sub_TR ' -FILE ' sub_TR ' one1.1D -tstim `cat tmp.1D` -numout ' n_tp ' > ' c '_' r '_' sd '.1D'...
                'rm -f tmp.1D'...
                '1deval -a Seed_Neur_' r '_' sd '.1D\ -b ' c '_' r '_' sd '.1D -expr "a*b" > Inter_neu' c '_' r '_' sd '.1D'...
                'waver -GAM -peak 1 -' TR ' ' sub_TR ' -input Inter_neu' c '_' r '_' sd '.1D -numout ' n_tp ' > Inter_hrf' c '_' r '_' sd '.1D'...
                '1dcat Inter_hrf' c '_' r '_' sd '.1D{0..$(2)} > Inter_ts' c '_' r '_' sd '.1D']);
      end
end


%% Catenate the two regressors across runs
unix(['cat Seed_ts?' sd 'D.1D > Seed_ts' sd '.1D']);
unix(['cat Inter_ts' c '?' sd '.1D > Inter_ts' c '_' sd '.1D']);

%% Re-run regression analysis by adding the two new regressors

3dDeconvolve -input pb04.$subj.r??.scale+tlrc.HEAD \
-polort A \
-mask full_mask.$subj+orig \
-num_stimts 13 \
-stim_times 1 stimuli/Allruns_stim_A_time.1D 'BLOCK(3,1)' \
-stim_label 1 A \
-stim_times 2 stimuli/Allruns_stim_B_time.1D 'BLOCK(3,1)' \
-stim_label 2 B \
-stim_times 3 stimuli/Allruns_stim_C_time.1D 'BLOCK(3,1)' \
-stim_label 3 C \
-stim_file 4 dfile.rall.1D'[0]' -stim_base 4 -stim_label 6 roll \
-stim_file 5 dfile.rall.1D'[1]' -stim_base 5 -stim_label 7 pitch \
-stim_file 6 dfile.rall.1D'[2]' -stim_base 6 -stim_label 8 yaw \
-stim_file 7 dfile.rall.1D'[3]' -stim_base 7 -stim_label 9 dS \
-stim_file 8 dfile.rall.1D'[4]' -stim_base 8 -stim_label 10 dL \
-stim_file 9 dfile.rall.1D'[5]' -stim_base 9 -stim_label 11 dP \
-stim_file 10 Seed_ts' sd '.1D -stim_label 10 Seed \
-stim_file 11 Inter_tsA' sd '.1D -stim_label 11 PPIA \
-stim_file 12 Inter_tsB' sd '.1D -stim_label 12 PPIB \
-stim_file 13 Inter_tsC' sd '.1D -stim_label 13 PPIC \
-rout -tout \
-bucket PPIstat' sd '