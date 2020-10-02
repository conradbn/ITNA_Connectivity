purge

%% --------- Digit LH, Dp --------- 
label = 'Digit_lh_Dp';
hemi = 'lh';
test = 1;
niters = 5000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/gPPI_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_gPPI_' label '.' hemi '.niml.dset'];

% Create catenated dateset based on statistic brik label
% [~,cmdout] = unix('3dinfo -label2index iDp#0_Coef *_proc/*.PPI/stats.PPI.*.Dp-Da.lh_REML.niml.dset');
brik1 = '19';
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.PPI/stats.PPI.*.Dp-Da.' hemi '_REML.niml.dset[' brik1 ']']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_gPPI_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

%% --------- Digit LH, Lp --------- 
label = 'Digit_lh_Lp';
hemi = 'lh';
test = 1;
niters = 5000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/gPPI_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_gPPI_' label '.' hemi '.niml.dset'];

% Create catenated dateset based on statistic brik label
% [~,cmdout] = unix('3dinfo -label2index iLa#0_Coef *_proc/*.PPI/stats.PPI.*.Dp-Da.lh_REML.niml.dset');
brik1 = '25';
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.PPI/stats.PPI.*.Dp-Da.' hemi '_REML.niml.dset[' brik1 ']']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_gPPI_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

%% --------- Digit LH, Dp - Da --------- 
label = 'Digit_lh_Dp_Da';
hemi = 'lh';
test = 2;
niters = 5000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/gPPI_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_gPPI_' label '.' hemi '.niml.dset'];

% Create catenated dateset based on statistic brik label
% [~,cmdout] = unix('3dinfo -label2index iDp-iDa#0_Coef *_proc/*.PPI/stats.PPI.*.Dp-Da.lh_REML.niml.dset');
% [~,cmdout] = unix('3dinfo -label2index iLp-iLa#0_Coef *_proc/*.PPI/stats.PPI.*.Dp-Da.lh_REML.niml.dset');
brik1 = '19';
brik2 = '22';
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.PPI/stats.PPI.*.Dp-Da.' hemi '_REML.niml.dset[' brik1 ']'...
     ' *_proc/*.PPI/stats.PPI.*.Dp-Da.' hemi '_REML.niml.dset[' brik2 ']']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_gPPI_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

%% --------- Digit LH, Dp - Lp --------- 
label = 'Digit_lh_Dp_Lp';
hemi = 'lh';
test = 2;
% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/gPPI_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_gPPI_' label '.' hemi '.niml.dset'];

% Create catenated dateset based on statistic brik label
% [~,cmdout] = unix('3dinfo -label2index iDp-iDa#0_Coef *_proc/*.PPI/stats.PPI.*.Dp-Da.lh_REML.niml.dset');
% [~,cmdout] = unix('3dinfo -label2index iLp-iLa#0_Coef *_proc/*.PPI/stats.PPI.*.Dp-Da.lh_REML.niml.dset');
brik1 = '19';
brik2 = '25';
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.PPI/stats.PPI.*.Dp-Da.' hemi '_REML.niml.dset[' brik1 ']'...
     ' *_proc/*.PPI/stats.PPI.*.Dp-Da.' hemi '_REML.niml.dset[' brik2 ']']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_gPPI_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)

%% --------- Digit LH, Dp-Da - Lp-La --------- 
label = 'Digit_lh_DpDa_LpLa';
hemi = 'lh';
test = 2;
niters = 5000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/gPPI_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_gPPI_' label '.' hemi '.niml.dset'];

% Create catenated dateset based on statistic brik label
% [~,cmdout] = unix('3dinfo -label2index iDp-iDa#0_Coef *_proc/*.PPI/stats.PPI.*.Dp-Da.lh_REML.niml.dset');
% [~,cmdout] = unix('3dinfo -label2index iLp-iLa#0_Coef *_proc/*.PPI/stats.PPI.*.Dp-Da.lh_REML.niml.dset');
brik1 = '34';
brik2 = '37';
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.PPI/stats.PPI.*.Dp-Da.' hemi '_REML.niml.dset[' brik1 ']'...
     ' *_proc/*.PPI/stats.PPI.*.Dp-Da.' hemi '_REML.niml.dset[' brik2 ']']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_gPPI_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)


%% --------- Digit LH, Dp-Da - Lp-La --------- 
label = 'Letter_lh_LpLa_DpDa';
hemi = 'lh';
test = 2;
niters = 5000;

% Create combined, covariate corrected residual dataset
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/')
out_dir = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/gPPI_Analysis/' label];
mkdir(out_dir);
out_dset = [out_dir '/allsubs_gPPI_' label '.' hemi '.niml.dset'];

% Create catenated dateset based on statistic brik label
% [~,cmdout] = unix('3dinfo -label2index iLp-iLa#0_Coef *_proc/*.PPI/stats.PPI.*.Lp-La.lh_REML.niml.dset');
% [~,cmdout] = unix('3dinfo -label2index iDp-iDa#0_Coef *_proc/*.PPI/stats.PPI.*.Lp-La.lh_REML.niml.dset');
brik1 = '37';
brik2 = '34';
unix(['3dTcat -overwrite -prefix ' out_dset ...
     ' *_proc/*.PPI/stats.PPI.*.Lp-La.' hemi '_REML.niml.dset[' brik1 ']'...
     ' *_proc/*.PPI/stats.PPI.*.Lp-La.' hemi '_REML.niml.dset[' brik2 ']']);

% Set filename for the output statistic dataset and number of iterations
out = [out_dir '/tfce_gPPI_' label '.' hemi '.' num2str(niters) 'iters.niml.dset'];

setup_run_TFCE(out_dset,niters,out,test,hemi)


%% Global function
function [] = setup_run_TFCE(out_dset,niters,out,test,hemi)

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
        
[~] = NFA_run_TFCE(surf_ds,vertices,faces,niters,out);
end

