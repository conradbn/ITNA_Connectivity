purge

%% --------- Digit LH, Lp --------- 
label = 'Digit_lh_Dp';
hemi = 'lh';
test = 1;
niters = 10000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

% Create catenated dateset
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zmap.Dp.niml.dset']);
% Create catenated dateset
unix(['3dMean -overwrite -prefix ' out_dset_mean ...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zmap.Dp.niml.dset']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

%% --------- Digit LH, Da --------- 
label = 'Digit_lh_Da';
hemi = 'lh';
test = 1;
niters = 10000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

% Create catenated dateset
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zmap.Da.niml.dset']);
% Create catenated dateset
unix(['3dMean -overwrite -prefix ' out_dset_mean ...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zmap.Da.niml.dset']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

%% --------- Digit LH, Dp-Da --------- 
label = 'Digit_lh_Dp-Da';
hemi = 'lh';
test = 1;
niters = 10000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

% Create catenated dateset
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zdiff.Dp-Da.niml.dset']);
% Create catenated dateset
unix(['3dMean -overwrite -prefix ' out_dset_mean ...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zdiff.Dp-Da.niml.dset']);

% Create catenated null dateset
out_dset_null = strrep(out_dset,'.niml.dset','_null.niml.dset');
unix(['3dTcat -overwrite -prefix ' out_dset_null ...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zdiff.Dp-Da_null.niml.dset']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi,out_dset_null)


%% --------- Digit RH, Dp --------- 
label = 'Digit_rh_Dp';
hemi = 'rh';
test = 1;
niters = 10000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

% Create catenated dateset
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.beta_series/*.Dp-Da.rh.beta_series_corr.' hemi '.Zmap.Dp.niml.dset']);
% Create catenated dateset
unix(['3dMean -overwrite -prefix ' out_dset_mean ...
     ' *_proc/*.beta_series/*.Dp-Da.rh.beta_series_corr.' hemi '.Zmap.Dp.niml.dset']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

%% --------- Digit RH, Da --------- 
label = 'Digit_rh_Da';
hemi = 'rh';
test = 1;
niters = 10000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

% Create catenated dateset
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.beta_series/*.Dp-Da.rh.beta_series_corr.' hemi '.Zmap.Da.niml.dset']);
% Create catenated dateset
unix(['3dMean -overwrite -prefix ' out_dset_mean ...
     ' *_proc/*.beta_series/*.Dp-Da.rh.beta_series_corr.' hemi '.Zmap.Da.niml.dset']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

%% --------- Digit RH, Dp-Da --------- 
label = 'Digit_rh_Dp-Da';
hemi = 'rh';
test = 1;
niters = 10000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

% Create catenated dateset
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.beta_series/*.Dp-Da.rh.beta_series_corr.' hemi '.Zdiff.Dp-Da.niml.dset']);
 % Create catenated dateset
unix(['3dMean -overwrite -prefix ' out_dset_mean ...
     ' *_proc/*.beta_series/*.Dp-Da.rh.beta_series_corr.' hemi '.Zdiff.Dp-Da.niml.dset']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)



                    %% --------- Letter LH, Lp --------- 
                    label = 'Letter_lh_Lp';
                    hemi = 'lh';
                    test = 1;
                    niters = 10000;

                    % Create combined, covariate corrected residual dataset
                    cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
                    out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
                    mkdir(out_dir);
                    out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
                    out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

                    % Create catenated dateset
                    unix(['3dTcat -overwrite -prefix ' out_dset ...
                         ' *_proc/*.beta_series/*.Lp-La.lh.beta_series_corr.' hemi '.Zmap.Lp.niml.dset']);
                    % Create catenated dateset
                    unix(['3dMean -overwrite -prefix ' out_dset_mean ...
                         ' *_proc/*.beta_series/*.Lp-La.lh.beta_series_corr.' hemi '.Zmap.Lp.niml.dset']);

                    % Set filename for the output statistic dataset and number of iterations
                    out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

                    setup_run_TFCE(out_dset,niters,out,test,hemi)

                    %% --------- Letter LH, La --------- 
                    label = 'Letter_lh_La';
                    hemi = 'lh';
                    test = 1;
                    niters = 10000;

                    % Create combined, covariate corrected residual dataset
                    cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
                    out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
                    mkdir(out_dir);
                    out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
                    out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

                    % Create catenated dateset
                    unix(['3dTcat -overwrite -prefix ' out_dset ...
                         ' *_proc/*.beta_series/*.Lp-La.lh.beta_series_corr.' hemi '.Zmap.La.niml.dset']);
                    % Create catenated dateset
                    unix(['3dMean -overwrite -prefix ' out_dset_mean ...
                         ' *_proc/*.beta_series/*.Lp-La.lh.beta_series_corr.' hemi '.Zmap.La.niml.dset']);

                    % Set filename for the output statistic dataset and number of iterations
                    out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

                    setup_run_TFCE(out_dset,niters,out,test,hemi)

                    %% --------- Letter LH, Lp-La --------- 
                    label = 'Letter_lh_Lp-La';
                    hemi = 'lh';
                    test = 1;
                    niters = 10000;

                    % Create combined, covariate corrected residual dataset
                    cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
                    out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
                    mkdir(out_dir);
                    out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
                    out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

                    % Create catenated dateset
                    unix(['3dTcat -overwrite -prefix ' out_dset ...
                         ' *_proc/*.beta_series/*.Lp-La.lh.beta_series_corr.' hemi '.Zdiff.Lp-La.niml.dset']);
                     % Create catenated dateset
                    unix(['3dMean -overwrite -prefix ' out_dset_mean ...
                         ' *_proc/*.beta_series/*.Lp-La.lh.beta_series_corr.' hemi '.Zdiff.Lp-La.niml.dset']);

                    % Set filename for the output statistic dataset and number of iterations
                    out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

                    setup_run_TFCE(out_dset,niters,out,test,hemi)


                    %% --------- Letter RH, Lp --------- 
                    label = 'Letter_rh_Lp';
                    hemi = 'rh';
                    test = 1;
                    niters = 10000;

                    % Create combined, covariate corrected residual dataset
                    cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
                    out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
                    mkdir(out_dir);
                    out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
                    out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

                    % Create catenated dateset
                    unix(['3dTcat -overwrite -prefix ' out_dset ...
                         ' *_proc/*.beta_series/*.Lp-La.rh.beta_series_corr.' hemi '.Zmap.Lp.niml.dset']);
                    % Create catenated dateset
                    unix(['3dMean -overwrite -prefix ' out_dset_mean ...
                         ' *_proc/*.beta_series/*.Lp-La.rh.beta_series_corr.' hemi '.Zmap.Lp.niml.dset']);

                    % Set filename for the output statistic dataset and number of iterations
                    out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

                    setup_run_TFCE(out_dset,niters,out,test,hemi)

                    %% --------- Letter RH, La --------- 
                    label = 'Letter_rh_La';
                    hemi = 'rh';
                    test = 1;
                    niters = 10000;

                    % Create combined, covariate corrected residual dataset
                    cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
                    out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
                    mkdir(out_dir);
                    out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
                    out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

                    % Create catenated dateset
                    unix(['3dTcat -overwrite -prefix ' out_dset ...
                         ' *_proc/*.beta_series/*.Lp-La.rh.beta_series_corr.' hemi '.Zmap.La.niml.dset']);
                    % Create catenated dateset
                    unix(['3dMean -overwrite -prefix ' out_dset_mean ...
                         ' *_proc/*.beta_series/*.Lp-La.rh.beta_series_corr.' hemi '.Zmap.La.niml.dset']);

                    % Set filename for the output statistic dataset and number of iterations
                    out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

                    setup_run_TFCE(out_dset,niters,out,test,hemi)

                    %% --------- Letter RH, Lp-La --------- 
                    label = 'Letter_rh_Lp-La';
                    hemi = 'rh';
                    test = 1;
                    niters = 10000;

                    % Create combined, covariate corrected residual dataset
                    cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
                    out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
                    mkdir(out_dir);
                    out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
                    out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

                    % Create catenated dateset
                    unix(['3dTcat -overwrite -prefix ' out_dset ...
                         ' *_proc/*.beta_series/*.Lp-La.rh.beta_series_corr.' hemi '.Zdiff.Lp-La.niml.dset']);
                     % Create catenated dateset
                    unix(['3dMean -overwrite -prefix ' out_dset_mean ...
                         ' *_proc/*.beta_series/*.Lp-La.rh.beta_series_corr.' hemi '.Zdiff.Lp-La.niml.dset']);

                    % Set filename for the output statistic dataset and number of iterations
                    out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

                    setup_run_TFCE(out_dset,niters,out,test,hemi)

        %% --------- Digit LH, Dp-Lp --------- 
        label = 'Digit_lh_Dp-Lp';
        hemi = 'lh';
        test = 1;
        niters = 1000;

        % Create combined, covariate corrected residual dataset
        cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
        out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
        mkdir(out_dir);
        out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
        out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

        % Create catenated dateset
        unix(['3dTcat -overwrite -prefix ' out_dset ...
             ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zdiff.Dp-Lp.niml.dset']);
         % Create catenated dateset
        unix(['3dMean -overwrite -prefix ' out_dset_mean ...
             ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zdiff.Dp-Lp.niml.dset']);

        % Set filename for the output statistic dataset and number of iterations
        out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

        setup_run_TFCE(out_dset,niters,out,test,hemi)
        
%% --------- Digit LH, Dp-Da (contrast) --------- 
label = 'Digit_lh_Dp-Da_paired';
hemi = 'lh';
test = 2;
niters = 10000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_BSC_' label '.' hemi '.niml.dset'];
out_dset_mean = [out_dir '/mean_BSC_' label '.' hemi '.niml.dset'];

% Create catenated dateset
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zmap.Dp.niml.dset'...
     ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zmap.Da.niml.dset']);
%  % Create catenated dateset
% unix(['3dMean -overwrite -prefix ' out_dset_mean ...
%      ' *_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' hemi '.Zdiff.Dp-Da.niml.dset']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_BSC_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

                    
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

