%% Run Threshold Free Cluster Enhancement (TFCE) significance testing for all pairwise contrasts
purge
data_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/AllSubs_dsets';
out_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Results';       

% Define the inputs for two-sample testing, specifying (N contasts x 4 cell array)
% 1) the label for the contrast
% 2) the group mask
% 3/4) the two datasets for comparison
%
% Each dataset contains a 2D matrix with every subject's connectivity map (surface vector) 
% catenated into one file.

input_strings = {
    % Structural Connectivity Contrasts
    'PairedTest_SC_DigLH_vs_LetLH_log+c','lh','ld141',...
    'GroupMask_DigLH_LetLH_ConsThr0.3_ld141.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset';
        
    'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf','lh','ld141',...
    'GroupMask_DigLH_DigRH_ConsThr0.3_on_LHsurf_ld141.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c_MAP2CON.niml.dset';
        
    'PairedTest_SC_DigLH_vs_DigRH_log+c_on_RHsurf','rh','ld141',...
    'GroupMask_DigLH_DigRH_ConsThr0.3_on_RHsurf_ld141.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c_MAP2CON.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c.niml.dset';
        
        
    % Functional Connectivity Contrasts - Digit vs Letter
    'PairedTest_FC_ALL_DigLH_vs_ALL_LetLH','lh','ld60',...
    'GroupMask_DigLH_LetLH_ld60.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.ALL.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.ALL.niml.dset';
        
        
    'PairedTest_FC_Dp_DigLH_vs_Da_DigLH','lh','ld60',...
    'LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Da.niml.dset';
        
    'PairedTest_FC_Lp_LetLH_vs_La_LetLH','lh','ld60',...
    'LitCoord_Letter_Pollack19_-42_-64_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.Lp.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.La.niml.dset';  
        
        
    'PairedTest_FC_Dp-Da_DigLH_vs_Lp-La_LetLH','lh','ld60',...
    'GroupMask_DigLH_LetLH_ld60.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zdiff.Lp-La.niml.dset';  
        
        
    % Functional Connectivity Contrasts - Task-level        
    'PairedTest_FC_DTask_DigLH_vs_LTask_DigLH','lh','ld60',...
    'LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.ALL_DTask.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.ALL_LTask.niml.dset';
            
    'PairedTest_FC_DTask_LetLH_vs_LTask_LetLH','lh','ld60',...
    'LitCoord_Letter_Pollack19_-42_-64_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.ALL_DTask.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.ALL_LTask.niml.dset';
        
        
    'PairedTest_FC_Dp_DigLH_vs_Lp_DigLH','lh','ld60',...
    'LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Lp.niml.dset';
        
    'PairedTest_FC_Dp_LetLH_vs_Lp_LetLH','lh','ld60',...
    'LitCoord_Letter_Pollack19_-42_-64_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.Dp.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.Lp.niml.dset';
        
        
    % Functional Connectivity Contrasts - Digit L vs Digit R
    'PairedTest_FC_ALL_DigLH_vs_DigRH_on_LHsurf','lh','ld60',...
    'GroupMask_DigLH_DigRH_on_LHsurf_ld60.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.ALL.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.ALL_MAP2CON.niml.dset';  
    
    'PairedTest_FC_ALL_DigLH_vs_DigRH_on_RHsurf','rh','ld60',...
    'GroupMask_DigLH_DigRH_on_LHsurf_ld60.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.ALL_MAP2CON.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.ALL.niml.dset';
    
    
    'PairedTest_FC_Dp_DigRH_vs_Da_DigRH_on_RHsurf','rh','ld60',...
    'LitCoord_Digit_Pollack19_54_-52_-14_std.60_rh.inflated.14mm_diam_INV.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Da.niml.dset';
        
        
    'PairedTest_FC_Dp-Da_DigLH_vs_DigRH_on_LHsurf','lh','ld60',...
    'GroupMask_DigLH_DigRH_on_LHsurf_ld60.niml.dset',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da_MAP2CON.niml.dset';       

    'PairedTest_FC_Dp-Da_DigLH_vs_DigRH_on_RHsurf','rh','ld60',...
    'GroupMask_DigLH_DigRH_on_RHsurf_ld60.niml.dset'...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da_MAP2CON.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da.niml.dset'};
   
   
%     'PairedTest_FC_Dp-Da_DigLH_vs_LetLH','lh','ld60',...
%         'GroupMask_DigLH_LetLH_ld60.niml.dset',...
%         'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset',...
%         'AllSubs_Lp-La.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset';
% 
%     'PairedTest_FC_Dp_DigLH_vs_Lp_LetLH','lh','ld60',...
%         'GroupMask_DigLH_LetLH_ld60.niml.dset',...
%         'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset',...
%         'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.Lp.niml.dset';
% 

%% Loop through and run each contrast through TFCE
cd(data_dir)
niters = 100000; % Total number of iterations for TFCE
test = 2; % All of the tests here are paired two-sample tests
for ii = 1:size(input_strings,1)
    % Get the label for this test 
    label = [out_dir '/' input_strings{ii,1}];
    hemi = input_strings{ii,2};
    density = input_strings{ii,3};
    mask_dset = input_strings{ii,4};
    
    % Create mean datasets for inspection/visualizaton
    unix(['3dTstat -mean -overwrite -prefix ' label '_SET1_MEAN '...
        input_strings{ii,5}]);
    unix(['3dTstat -mean -overwrite -prefix ' label '_SET2_MEAN '...
        input_strings{ii,6}]);
    
    % Combine all data into one input file
    unix(['3dTcat -overwrite -prefix ' label '_ALLDATA '...
        input_strings{ii,5} ' ' input_strings{ii,6}]);
    
    in_dset = [label '_ALLDATA.niml.dset'];
    
    % Set filename for the output statistic dataset and number of iterations
    out = [label '_TFCE_Zscore_' num2str(niters) 'iters_MASK.niml.dset'];
    
    % Specify no input null (i.e. use the program's default method)
    out_dset_null = '';
    
    % Call global function
    setup_run_TFCE(in_dset,mask_dset,niters,out,test,hemi,density,out_dset_null)
end


%% Run conjunction analyses
% This takes all the paired contrasts, modifies their inputs, then runs the
% conjunction on their thresholded overlap maps
cd(data_dir)
niters = 100000; % Total number of iterations for TFCE
test = 1; % All of the tests here are one-sample tests

for ii = 1:3%size(input_strings,1)
    % Get the label for this test 
    label = [out_dir '/' input_strings{ii,1}];
    hemi = input_strings{ii,2};
    density = input_strings{ii,3};
    mask_dset = input_strings{ii,4};
    
    % Change label to indicate conjunction
    label = strrep(label,'vs','CONJ');
    
    % Get thresholded data for input into conjunction TFCE (1-sample test
    % on the overlap maps)
    if contains(label,'PairedTest_SC')
        % Threhsold stuctural connectivity maps
        cutoff = -6;
        thr_name1 = strrep(input_strings{ii,5},'.niml.dset','.thr-6.niml.dset');
        ds = afni_niml_readsimple(input_strings{ii,5});
        ds.data(ds.data > cutoff) = 1;
        ds.data(ds.data <= cutoff) = 0;
        afni_niml_writesimple(ds,thr_name1);
        
        thr_name2 = strrep(input_strings{ii,6},'.niml.dset','.thr-6.niml.dset');
        ds = afni_niml_readsimple(input_strings{ii,6});
        ds.data(ds.data > cutoff) = 1;
        ds.data(ds.data <= cutoff) = 0;
        afni_niml_writesimple(ds,thr_name2);
        
        % Also use the conjunction mask (overlap from consistency
        % thresholding)
        mask_dset = strrep(mask_dset,'DigLH','DigLH_CONJ');
        
    elseif contains(label,'PairedTest_FC') && contains(input_strings{ii,5},'Zmap')
        thr_name1 = strrep(input_strings{ii,5},'Zmap','Pval_thr005');
        thr_name2 = strrep(input_strings{ii,6},'Zmap','Pval_thr005');
    elseif contains(label,'PairedTest_FC') && contains(input_strings{ii,5},'Zdiff')
        thr_name1 = strrep(input_strings{ii,5},'Zdiff','Zdiff_Pval_thr01');
        thr_name2 = strrep(input_strings{ii,6},'Zdiff','Zdiff_Pval_thr01');
    end
    

    % Create mean datasets for inspection/visualizaton
    unix(['3dTstat -mean -overwrite -prefix ' label '_SET1_MEAN '...
        thr_name1]);
    unix(['3dTstat -mean -overwrite -prefix ' label '_SET2_MEAN '...
        thr_name2]);
    
    % Get subject-level overlap maps
    unix(['3dcalc -overwrite -a ' thr_name1 ' -b ' thr_name2...
         ' -prefix ' label '_ALLDATA -expr "and(a,b)"']); 
        
    in_dset = [label '_ALLDATA.niml.dset'];
    
    % Set filename for the output statistic dataset and number of iterations
    out = [label '_TFCE_Zscore_' num2str(niters) 'iters_MASK.niml.dset'];
    
    % Specify no input null (i.e. use the program's default method)
    out_dset_null = '';
    
    
    % Create count map (to potentially use in place of TFCE result)
    ds = afni_niml_readsimple(in_dset);
    ds.data = sum(ds.data,2)./size(ds.data,2);
    ds_mask = afni_niml_readsimple(mask_dset);
    ds.data = ds.data .* ds_mask.data;
    out_count = [label '_ALLDATA_count.niml.dset'];
    afni_niml_writesimple(ds,out_count);
    
    % Call global function
    %setup_run_TFCE(in_dset,mask_dset,niters,out,test,hemi,density,out_dset_null)
    
end
    
    

%% Global function
function [] = setup_run_TFCE(in_dset,mask_dset,niters,out,test,hemi,density,out_dset_null)

% Load in the subject data
surf_ds = afni_niml_readsimple(in_dset);
% Mask the subject data
mask_ds = afni_niml_readsimple(['/Users/nbl_imac2/Documents/GitHub/ITNA_Connectivity/roi_creation/rois/' mask_dset]);
surf_ds.data = surf_ds.data .* mask_ds.data;

% Convert to COSMOMVPA formatted structure
surf_ds.labels = num2cell(1:size(surf_ds.data,2));
surf_ds = cosmo_surface_dataset(surf_ds);

% Get faces and vertices info from MNI surface gii file
if strcmp(density,'ld60')
    surf_gii = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/suma_MNI152_2009/std.60.' hemi '.smoothwm.gii'];
elseif strcmp(density,'ld141')
    surf_gii = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/suma_MNI152_2009/std.141.' hemi '.smoothwm.gii'];
end

[vertices,faces]=surfing_read(surf_gii);
faces = double(faces);
vertices = double(vertices);

% Set up surface dataset structure, to prepare it for Cosmo TFCE
surf_ds.fa.center_ids = surf_ds.fa.node_indices;
if test == 2
    % Two-sample test
    nsubs = size(surf_ds.samples,1)/2;
    surf_ds.sa.chunks = [(1:nsubs)';(1:nsubs)'];
    surf_ds.sa.targets = [ones(nsubs,1);2*ones(nsubs,1)];
elseif test == 1
    % One-sample test (vs 0)
    nsubs = size(surf_ds.samples,1);
    surf_ds.sa.chunks = (1:nsubs)';
    surf_ds.sa.targets = ones(nsubs,1);
end

% Set up null surface data structure
surf_dset_null = cell(1,100);
if strcmp(out_dset_null,'') ~= 1
    null = afni_niml_readsimple(out_dset_null);
    ds = surf_ds;
    for ii = 1:100
        ds.samples = null.data(:,ii:100:size(null.data,2))';
        surf_dset_null(ii) = {ds};
    end
else
    surf_dset_null = [];
end

% Run TFCE
[~] = ITNA_run_TFCE(surf_ds,vertices,faces,niters,out,surf_dset_null);

% Also get the uncorrected t-stats and z-scores for verification 
% (NOTE this is adapted from last part of ITNA_run_TFCE)
stat_ds = cosmo_stat(surf_ds,'t','z');
tstat_ds = cosmo_stat(surf_ds,'t');
fprintf('TFCE z-score dataset\n');
cosmo_disp(stat_ds);
nfeatures=size(stat_ds.samples,2);
percentiles=(1:nfeatures)/nfeatures*100;
plot(percentiles,sort(stat_ds.samples))
title('sorted TFCE z-scores');
xlabel('feature percentile');
ylabel('z-score');
nvertices=size(vertices,1);
disp_opt=struct();
disp_opt.DataRange=[-2 2];
DispIVSurf(vertices,faces,1:nvertices,stat_ds.samples',0,disp_opt);

% Write uncorrected stats to surface files
out_unc = strrep(out,'.niml.dset','_uncorrZscr.niml.dset');
cosmo_map2surface(stat_ds,out_unc);
out_uncT = strrep(out,'.niml.dset','_uncorrT.niml.dset');
cosmo_map2surface(tstat_ds,out_uncT);

% Get Bayes Factors and write to surface file
X = surf_ds.samples(1:nsubs,:);
Y = surf_ds.samples(nsubs+1:end,:);
bf10_ds = tstat_ds;
parfor nn = 1:size(X,2)
    x = X(:,nn)';
    y = Y(:,nn)';
    bf_xy(nn) = bf.ttest(x,y);   
end
bf10_ds.samples(1,:) = bf_xy;
out_bf10 = strrep(out,'.niml.dset','_uncorrBF10.niml.dset');
cosmo_map2surface(bf10_ds,out_bf10);

end
