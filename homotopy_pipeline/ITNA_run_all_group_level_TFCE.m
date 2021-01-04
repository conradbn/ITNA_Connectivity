%% Run Threshold Free Cluster Enhancement (TFCE) significance testing for all pairwise contrasts
purge
niters = 100000; % Total number of iterations for TFCE
test = 2; % All of the tests here are paired two-sample tests
data_dir = '/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/AllSubs_dsets';
out_dir = '/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/Results';       

% Define the inputs for two-sample testing, specifying the label for the
% contrast then the two datasets for comparison (N contasts x 3 cell
% array). Each data file is a 4D dataset with every subject's connectivity
% map catenated into one file.

input_strings = {
    'PairedTest_SC_DigLH_vs_LetLH_log+c','lh','ld141',...
        'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset';
        
    'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf','lh','ld141',...
        'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c_MAP2CON.niml.dset';
        
    'PairedTest_SC_DigLH_vs_DigRH_log+c_on_RHsurf','rh','ld141',...
        'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log+c_MAP2CON.niml.dset',...
        'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c.niml.dset';
        
    'PairedTest_FC_Dp_DigLH_vs_LetLH','lh','ld60',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zmap.Dp.niml.dset';
        
    'PairedTest_FC_Dp-Da_DigLH_vs_LetLH','lh','ld60',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset',...
        'AllSubs_Lp-La.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset';
        
    'PairedTest_FC_Dp_DigLH_vs_DigRH_onLHsurf','lh','ld60',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp_MAP2CON.niml.dset';
        
    'PairedTest_FC_Dp-Da_DigLH_vs_DigRH_on_LHsurf','lh','ld60',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da_MAP2CON.niml.dset';
    
    'PairedTest_FC_Dp_DigLH_vs_DigRH_on_RHsurf','rh','ld60',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp_MAP2CON.niml.dset';
        
    'PairedTest_FC_Dp-Da_DigLH_vs_DigRH_on_RHsurf','rh','ld60',...
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da_MAP2CON.niml.dset',...
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da.niml.dset'};
    
% Linear models where connectivity is related to selectivity    
%     'SC_DigitLH_x_DigitSel'...
%     'FC_DigitLH_x_DigitSel'...  


% Loop through and run each contrast through TFCE
cd(data_dir)
for ii = 1:size(input_strings,1)
    % Get the label for this test 
    label = [out_dir '/' input_strings{ii,1}];
    hemi = input_strings{ii,2};
    density = input_strings{ii,3};
    
    % Create mean datasets for inspection/visualizaton
    unix(['3dTstat -mean -overwrite -prefix ' label '_SET1_MEAN '...
        input_strings{ii,4}]);
    unix(['3dTstat -mean -overwrite -prefix ' label '_SET2_MEAN '...
        input_strings{ii,5}]);
    
    % Combine all data into one input file
    unix(['3dTcat -overwrite -prefix ' label '_ALLDATA '...
        input_strings{ii,4} ' ' input_strings{ii,5}]);
    
    in_dset = [label '_ALLDATA.niml.dset'];
    
    % Set filename for the output statistic dataset and number of iterations
    out = [label '_TFCE_Zscore_' num2str(niters) 'iters.niml.dset'];
    
    % Specify no input null (i.e. use the program's default method)
    out_dset_null = '';
    
    % Call global function
    setup_run_TFCE(in_dset,niters,out,test,hemi,density,out_dset_null)
end


%% Global function
function [] = setup_run_TFCE(in_dset,niters,out,test,hemi,density,out_dset_null)

% Load in the subject data
surf_ds = afni_niml_readsimple(in_dset);
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

end
