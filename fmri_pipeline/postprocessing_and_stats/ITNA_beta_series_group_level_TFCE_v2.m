%% Run Threshold free cluster enhancement significance testing
purge
niters = 100000;
test = 1;

% Set input strings for each run
%                 seed   seed_hemi  data_hemi  condition    dtype
input_strings = {'Dp-Da'   'lh'       'lh'      'Dp-Da'    'Zdiff';...
                 'Dp-Da'   'lh'       'lh'      'Dp'       'Zmap' ;...
                 'Dp-Da'   'lh'       'lh'      'Da'       'Zmap' ;...
                 'Dp-Da'   'lh'       'lh'      'Dp-Lp'    'Zdiff';...
                 'Dp-Da'   'lh'       'lh'      'Lp'       'Zmap' ;...
                 'Dp-Da'   'lh'       'lh'      'La'       'Zmap' ;...
                 'Lp-La'   'lh'       'lh'      'Lp-La'    'Zdiff';...
                 'Lp-La'   'lh'       'lh'      'Dp'       'Zmap' ;...
                 'Lp-La'   'lh'       'lh'      'Da'       'Zmap' ;...
                 'Lp-La'   'lh'       'lh'      'Dp-Lp'    'Zdiff';...
                 'Lp-La'   'lh'       'lh'      'Lp'       'Zmap' ;...
                 'Lp-La'   'lh'       'lh'      'La'       'Zmap' ;...
                 'Dp-Da'   'rh'       'rh'      'Dp-Da'    'Zdiff';...
                 'Dp-Da'   'rh'       'rh'      'Dp'       'Zmap' ;...
                 'Dp-Da'   'rh'       'rh'      'Da'       'Zmap' ;...
                 'Dp-Da'   'rh'       'rh'      'Dp-Lp'    'Zdiff';...
                 'Dp-Da'   'rh'       'rh'      'Lp'       'Zmap' ;...
                 'Dp-Da'   'rh'       'rh'      'La'       'Zmap' ;...
                 'Lp-La'   'rh'       'rh'      'Lp-La'    'Zdiff';...
                 'Lp-La'   'rh'       'rh'      'Dp'       'Zmap' ;...
                 'Lp-La'   'rh'       'rh'      'Da'       'Zmap' ;...
                 'Lp-La'   'rh'       'rh'      'Dp-Lp'    'Zdiff';...
                 'Lp-La'   'rh'       'rh'      'Lp'       'Zmap' ;...
                 'Lp-La'   'rh'       'rh'      'La'       'Zmap' ;...
                 'Dp-Da'   'lh'       'lh'      'Lp-La'    'Zdiff'};
                 

% Loop through each group of settings             
for ii = 1:size(input_strings,1)
    % Set strings for this run
    seed = input_strings{ii,1};
    hemi_seed = input_strings{ii,2};
    hemi_data = input_strings{ii,3};
    condition = input_strings{ii,4};
    dtype = input_strings{ii,5};
    label = [seed '_' hemi_seed '_seed.' condition '_' hemi_data '_data'];
    
    % Set up data paths
    cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
    out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
    mkdir(out_dir);
    out_dset = [out_dir '/allsubs_BSC_' label '.' hemi_data '_data.niml.dset'];
    out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi_data '_data.niml.dset'];
       
    % Create catenated dateset
    unix(['3dTcat -overwrite -prefix ' out_dset ...
        ' *_proc/*.beta_series/*.' seed '.' hemi_seed '.beta_series_corr.' hemi_data '.' dtype '.' condition '.niml.dset']);
    % Create mean dateset
    unix(['3dMean -overwrite -prefix ' out_dset_mean ...
        ' *_proc/*.beta_series/*.' seed '.' hemi_seed '.beta_series_corr.' hemi_data '.' dtype '.' condition '.niml.dset']);
    
    % Create catenated null dateset
    out_dset_null = '';
%     out_dset_null = strrep(out_dset,'.niml.dset','_null.niml.dset');
%     unix(['3dTcat -overwrite -prefix ' out_dset_null ...
%         ' *_proc/*.beta_series/*.' seed '.' hemi_seed '.beta_series_corr.' hemi_data '.' dtype '.' condition '_null.niml.dset']);
    
    % Set filename for the output statistic dataset and number of iterations
    out = [out_dir '/tfce_BSC_' label '_seed.' hemi_data '_data.' num2str(niters) 'iters.niml.dset'];
%     out = [out_dir '/tfce_BSC_' label '_seed.' hemi_data '_data.' num2str(niters) 'iters.with_NULL.niml.dset'];
    
    % Call global function
    setup_run_TFCE(out_dset,niters,out,test,hemi_data,out_dset_null)
end


%% Global function
function [] = setup_run_TFCE(out_dset,niters,out,test,hemi,out_dset_null)

% Load in the subject data
surf_ds = afni_niml_readsimple(out_dset);
surf_ds.labels = num2cell(1:size(surf_ds.data,2));
surf_ds = cosmo_surface_dataset(surf_ds);

% Get faces and vertices info from MNI surface gii file
[vertices,faces]=surfing_read(['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/suma_MNI152_2009/std.60.' hemi '.smoothwm.gii']);
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
