%% Run Threshold free cluster enhancement significance testing
purge
niters = 100000;
test = 1;

% Set input strings for each run
%                 seed   seed_hemi  data_hemi  condition    dtype
input_strings = {'Dp-Da'   'lh'       'lh'      'Dp'       'Pval_thr005' ;...
                 'Lp-La'   'lh'       'lh'      'Dp'       'Pval_thr005' ;...
                 'Dp-Da'   'lh'       'lh'      'Da'       'Pval_thr005' ;...
                 'Lp-La'   'lh'       'lh'      'Da'       'Pval_thr005' ;...
                 'Dp-Da'   'lh'       'lh'      'Lp'       'Pval_thr005' ;...
                 'Lp-La'   'lh'       'lh'      'Lp'       'Pval_thr005' ;...
                 'Dp-Da'   'lh'       'lh'      'La'       'Pval_thr005' ;...
                 'Lp-La'   'lh'       'lh'      'La'       'Pval_thr005' ;...
                 'Dp-Da'   'rh'       'rh'      'Dp'       'Pval_thr005' ;...
                 'Lp-La'   'rh'       'rh'      'Dp'       'Pval_thr005' ;...
                 'Dp-Da'   'rh'       'rh'      'Da'       'Pval_thr005' ;...
                 'Lp-La'   'rh'       'rh'      'Da'       'Pval_thr005' ;...
                 'Dp-Da'   'rh'       'rh'      'Lp'       'Pval_thr005' ;...
                 'Lp-La'   'rh'       'rh'      'Lp'       'Pval_thr005' ;...
                 'Dp-Da'   'rh'       'rh'      'La'       'Pval_thr005' ;...
                 'Lp-La'   'rh'       'rh'      'La'       'Pval_thr005' };
                 
                 
             
% Loop through each group of settings             
for ii =1:2:size(input_strings,1)
    i1 = ii;
    i2 = ii+1;
    % Set strings for this run
    seed1 = input_strings{i1,1};
    hemi_seed1 = input_strings{i1,2};
    hemi_data1 = input_strings{i1,3};
    condition1 = input_strings{i1,4};
    dtype1 = input_strings{i1,5};
    label1 = [seed1 '_' hemi_seed1 '_seed.' condition1 '_' hemi_data1 '_data'];
    
    % Set strings for this run
    seed2 = input_strings{i2,1};
    hemi_seed2 = input_strings{i2,2};
    hemi_data2 = input_strings{i2,3};
    condition2 = input_strings{i2,4};
    dtype2 = input_strings{i2,5};
    label2 = [seed2 '_' hemi_seed2 '_seed.' condition2 '_' hemi_data2 '_data'];
    
    % Set up data paths
    cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
    out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label1 '.CONJ.' label2];
    mkdir(out_dir);
    
    % Create catenated datesets
    out_dset = [out_dir '/allsubs_BSC_' label1 '.BIN.niml.dset'];
    unix(['3dTcat -overwrite -prefix ' out_dset ...
        ' *_proc/*.beta_series/*.' seed1 '.' hemi_seed1 '.beta_series_corr.' hemi_data1 '.' dtype1 '.' condition1 '.niml.dset']);
    out_dset = [out_dir '/allsubs_BSC_' label2 '.BIN.niml.dset'];
    unix(['3dTcat -overwrite -prefix ' out_dset ...
        ' *_proc/*.beta_series/*.' seed2 '.' hemi_seed2 '.beta_series_corr.' hemi_data2 '.' dtype2 '.' condition2 '.niml.dset']);
       
    % Create mean datesets
    out_dset = [out_dir '/mean_BSC_' label1 '.BIN.niml.dset'];
    unix(['3dMean -overwrite -prefix ' out_dset ...
        ' *_proc/*.beta_series/*.' seed1 '.' hemi_seed1 '.beta_series_corr.' hemi_data1 '.' dtype1 '.' condition1 '.niml.dset']);
    out_dset = [out_dir '/mean_BSC_' label2 '.BIN.niml.dset'];
    unix(['3dMean -overwrite -prefix ' out_dset ...
        ' *_proc/*.beta_series/*.' seed2 '.' hemi_seed2 '.beta_series_corr.' hemi_data2 '.' dtype2 '.' condition2 '.niml.dset']);
    
    % Create overlap (conjunction) dataset
    dset1 = afni_niml_readsimple([out_dir '/allsubs_BSC_' label1 '.BIN.niml.dset']);
    dset2 = afni_niml_readsimple([out_dir '/allsubs_BSC_' label2 '.BIN.niml.dset']);
    dset_union = dset1;
    dset_union.data = dset1.data == 1 & dset2.data == 1;
    out_dset_union = [out_dir '/allsubs_BSC_' label1 '.CONJ.' label2 '.niml.dset'];
    afni_niml_writesimple(dset_union,out_dset_union); 
    
    % Create catenated null dateset
    out_dset_null = '';
%     out_dset_null = strrep(out_dset,'.niml.dset','_null.niml.dset');
%     unix(['3dTcat -overwrite -prefix ' out_dset_null ...
%         ' *_proc/*.beta_series/*.' seed '.' hemi_seed '.beta_series_corr.' hemi_data '.' dtype '.' condition '_null.niml.dset']);
    
    % Set filename for the output statistic dataset and number of iterations
    out = [out_dir '/tfce_BSC_' label1 '.CONJ.' label2 '.' num2str(niters) '.niml.dset'];
%     out = [out_dir '/tfce_BSC_' label '_seed.' hemi_data '_data.' num2str(niters) 'iters.with_NULL.niml.dset'];
    
    % Call global function
    setup_run_TFCE(out_dset_union,niters,out,test,hemi_data1,out_dset_null)
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
