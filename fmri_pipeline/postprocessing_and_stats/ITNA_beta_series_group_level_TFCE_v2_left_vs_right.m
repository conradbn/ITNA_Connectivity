%% Run Threshold free cluster enhancement significance testing
purge
niters = 100000;
test = 2;

% % Set input strings for each run
% input_strings = {'Dp-Da'   'lh'       'lh'      'Dp-Da'    'Zdiff';... % 1 
%                  'Lp-La'   'lh'       'lh'      'Lp-La'    'Zdiff';... % 2
%                  'Dp-Da'   'lh'       'lh'      'Dp'       'Zmap' ;... % 1
%                  'Lp-La'   'lh'       'lh'      'Dp'       'Zmap' ;... % 2
%                  'Dp-Da'   'lh'       'lh'      'Da'       'Zmap' ;... % 1
%                  'Lp-La'   'lh'       'lh'      'Da'       'Zmap' ;... % 2
%                  'Dp-Da'   'lh'       'lh'      'Lp'       'Zmap' ;... % 1
%                  'Lp-La'   'lh'       'lh'      'Lp'       'Zmap' ;... % 2
%                  'Dp-Da'   'lh'       'lh'      'La'       'Zmap' ;... % 1
%                  'Lp-La'   'lh'       'lh'      'La'       'Zmap' ;... % 2
%                  'Dp-Da'   'rh'       'rh'      'Dp-Da'    'Zdiff';... % 1 
%                  'Lp-La'   'rh'       'rh'      'Lp-La'    'Zdiff';... % 2
%                  'Dp-Da'   'rh'       'rh'      'Dp'       'Zmap' ;... % 1
%                  'Lp-La'   'rh'       'rh'      'Dp'       'Zmap' ;... % 2
%                  'Dp-Da'   'rh'       'rh'      'Da'       'Zmap' ;... % 1
%                  'Lp-La'   'rh'       'rh'      'Da'       'Zmap' ;... % 2
%                  'Dp-Da'   'rh'       'rh'      'Lp'       'Zmap' ;... % 1
%                  'Lp-La'   'rh'       'rh'      'Lp'       'Zmap' ;... % 2
%                  'Dp-Da'   'rh'       'rh'      'La'       'Zmap' ;... % 1
%                  'Lp-La'   'rh'       'rh'      'La'       'Zmap' ;... % 2
%                  
%                  };
% Loop through each group of settings             
% for ii = 1:2:size(input_strings,1)
%     i1 = ii;
%     i2 = ii+1;
%     % Set strings for this run
%     seed1 = input_strings{i1,1};
%     hemi_seed1 = input_strings{i1,2};
%     hemi_data1 = input_strings{i1,3};
%     condition1 = input_strings{i1,4};
%     dtype1 = input_strings{i1,5};
%     label1 = [seed1 '_' hemi_seed1 '_seed.' condition1 '_' hemi_data1 '_data'];
%     
%     % Set strings for this run
%     seed2 = input_strings{i2,1};
%     hemi_seed2 = input_strings{i2,2};
%     hemi_data2 = input_strings{i2,3};
%     condition2 = input_strings{i2,4};
%     dtype2 = input_strings{i2,5};
%     label2 = [seed2 '_' hemi_seed2 '_seed.' condition2 '_' hemi_data2 '_data'];
%     
%     % Set up data paths
%     cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
%     out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label1 '_MINUS_' label2];
%     mkdir(out_dir);
%     % Create catenated dateset
%     out_dset = [out_dir '/allsubs_BSC_' label1 '.MINUS.' label2 '.niml.dset'];
%     unix(['3dTcat -overwrite -prefix ' out_dset ...
%         ' *_proc/*.beta_series/*.' seed1 '.' hemi_seed1 '.beta_series_corr.' hemi_data1 '.' dtype1 '.' condition1 '.niml.dset'...
%         ' *_proc/*.beta_series/*.' seed2 '.' hemi_seed2 '.beta_series_corr.' hemi_data2 '.' dtype2 '.' condition2 '.niml.dset']);
%     % Create mean dateset
%     out_dset_mean = [out_dir '/mean_BSC_' label1 '.MINUS.' label2 '.niml.dset'];
%     unix(['3dMean -overwrite -prefix ' out_dset_mean ...
%         ' *_proc/*.beta_series/*.' seed '.' hemi_seed '.beta_series_corr.' hemi_data '.' dtype '.' condition '.niml.dset']);
    
%     % Create catenated null dateset
%     out_dset_null = strrep(out_dset,'.niml.dset','_null.niml.dset');
%     unix(['3dTcat -overwrite -prefix ' out_dset_null ...
%         ' *_proc/*.beta_series/*.' seed1 '.' hemi_seed1 '.beta_series_corr.' hemi_data1 '.' dtype1 '.' condition1 '_null.niml.dset'...
%         ' *_proc/*.beta_series/*.' seed2 '.' hemi_seed2 '.beta_series_corr.' hemi_data2 '.' dtype2 '.' condition2 '_null.niml.dset']);
%     
%     % Set filename for the output statistic dataset and number of iterations
%     out = [out_dir '/tfce_BSC_' label1 '.MINUS.' label2 '.' num2str(niters) 'iters.with_NULL.niml.dset'];
%     
%     out_dset_null = '';
%     % Call global function
%     setup_run_TFCE(out_dset,niters,out,test,hemi_data1,out_dset_null)
% end



purge
niters = 1000;
test = 2; 
cd('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/AllSubs_dsets');
input_strings = {'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset'
                 'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp_MAPPED2CONTRA.niml.dset'
                 'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset'
                 'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da_MAPPED2CONTRA.niml.dset'};
             
out_dir = '/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/Results';             
out_dset = [out_dir '/AllSubs_Dp_Da.lhtolh_vs_rhtorh_Zmap.Dp.niml.dset'];
unix(['3dTcat -overwrite -prefix ' out_dset ...
      ' ' input_strings{1} ' ' input_strings{2}]);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_BSC_' num2str(niters) 'iters.Dp_Da.lhtolh_vs_rhtorh_Zmap.Dp.niml.dset'];
out_dset_null = '';
hemi_data1 = 'lh';

% Call global function
setup_run_TFCE(out_dset,niters,out,test,hemi_data1,out_dset_null)


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
