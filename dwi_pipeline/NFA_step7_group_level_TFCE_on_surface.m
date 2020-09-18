purge

%% --------- LEFT HEMISPHERE --------- 
% Create combined, covariate corrected residual dataset
cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/TFCE')
out = 'allsubs_lh_Digit_Letter.lh.TDI_ends.norm.6mm.log+c.niml.dset';
unix(['3dTcat -overwrite -prefix ' out ...
     ' ../*/*/resid_from_cov_Digit_TDI_ends_surf_lh_6mm_log+c_3dttest++.lh.mask.niml.dset'...
     ' ../*/*/resid_from_cov_Letter_TDI_ends_surf_lh_6mm_log+c_3dttest++.lh.mask.niml.dset']);

% Load in the subject data
surf_ds = afni_niml_readsimple(out);
surf_ds.labels = num2cell(1:size(surf_ds.data,2));
surf_ds = cosmo_surface_dataset(surf_ds);

% Get faces and vertices info from MNI surface gii file
[vertices,faces]=surfing_read('../suma_MNI152_2009/std.141.lh.smoothwm.gii');
faces = double(faces);
vertices = double(vertices);

% Set up surface dataset structure, to prepare it for Cosmo TFCE
nsubs = size(surf_ds.samples,1)/2;
surf_ds.fa.center_ids = surf_ds.fa.node_indices;
surf_ds.sa.chunks = [(1:nsubs)';(1:nsubs)'];
surf_ds.sa.targets = [ones(nsubs,1);2*ones(nsubs,1)];
        % surf_ds.sa.chunks = (1:nsubs)';
        % surf_ds.sa.targets = ones(nsubs,1);

% Set filename for the output statistic dataset and number of iterations
out = 'tfce_lh_D-L_log+c.lh.50k.niml.dset';
niters = 50000;

% Mask the dataset to reduce to only areas of "consistent" connectivity to
% digit OR letter ROIs
mask_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Data';
m1 = [mask_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'];
m2 = [mask_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'];
m3 = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Masks/lh_rois.niml.dset'];

m1s = afni_niml_readsimple(m1);
m2s = afni_niml_readsimple(m2);
m3s = afni_niml_readsimple(m3);
ind_keep = m1s.data == 1 | m2s.data == 1;
ind_keep(m3s.data(:,2)==1) = 0;
surf_ds.samples(:,~ind_keep) = NaN;

surf_ds_lh = surf_ds;
vertices_lh = vertices;
faces_lh = faces;
niters_lh = niters;
out_lh = out;

clearvars -except surf_ds_lh vertices_lh faces_lh niters_lh out_lh


%% --------- RIGHT HEMISPHERE --------- 
% Create combined, covariate corrected residual dataset
cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/TFCE')
out = 'allsubs_rh_Digit_Letter.rh.TDI_ends.norm.6mm_log+c.niml.dset';
unix(['3dTcat -overwrite -prefix ' out ...
     ' ../*/*/resid_from_cov_Digit_TDI_ends_surf_rh_6mm_log+c_3dttest++.rh.mask.niml.dset'...
     ' ../*/*/resid_from_cov_Letter_TDI_ends_surf_rh_6mm_log+c_3dttest++.rh.mask.niml.dset']);
 
% Load in the subject data
surf_ds = afni_niml_readsimple(out);
surf_ds.labels = num2cell(1:size(surf_ds.data,2));
surf_ds = cosmo_surface_dataset(surf_ds);
 
% Get faces and vertices info from MNI surface gii file
[vertices,faces]=surfing_read('../suma_MNI152_2009/std.141.rh.smoothwm.gii');
faces = double(faces);
vertices = double(vertices);

% Set up surface dataset structure, to prepare it for Cosmo TFCE
nsubs = size(surf_ds.samples,1)/2;
surf_ds.fa.center_ids = surf_ds.fa.node_indices;
surf_ds.sa.chunks = [(1:nsubs)';(1:nsubs)'];
surf_ds.sa.targets = [ones(nsubs,1);2*ones(nsubs,1)];
        % surf_ds.sa.chunks = (1:nsubs)';
        % surf_ds.sa.targets = ones(nsubs,1);

% Set filename for the output statistic dataset and number of iterations
out = 'tfce_rh_D-L_log+c.rh.50k.niml.dset';
niters = 50000;

% Mask the dataset to reduce to only areas of "consistent" connectivity to
% digit OR letter ROIs
mask_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Data';
m1 = [mask_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_Dp-Da.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'];
m2 = [mask_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_Lp-La.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'];
m3 = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Masks/rh_rois.niml.dset';

m1s = afni_niml_readsimple(m1);
m2s = afni_niml_readsimple(m2);
m3s = afni_niml_readsimple(m3);
ind_keep = m1s.data == 1 | m2s.data == 1;
ind_keep(m3s.data(:,2)==1) = 0;
surf_ds.samples(:,~ind_keep) = NaN;


surf_ds_rh = surf_ds;
vertices_rh = vertices;
faces_rh = faces;
niters_rh = niters;
out_rh = out;

%% Loop and run TFCE 
parfor ii = 1:2
    if ii ==1
        [~] = NFA_run_TFCE(surf_ds_lh,vertices_lh,faces_lh,niters_lh,out_lh);
    elseif ii ==2
        [~] = NFA_run_TFCE(surf_ds_rh,vertices_rh,faces_rh,niters_rh,out_rh);
    end
end

%% Create masks of significant difference
cutoff = 1.96;
out = out_lh;
h = afni_niml_readsimple(out_lh);
h.data(abs(h.data) <= cutoff) = 0;
h.data(h.data ~= 0) = 1;
afni_niml_writesimple(h,['mask_p05_' out]);
diff = surf_ds_lh.samples(1:29,:) - surf_ds_lh.samples(30:58,:);
diff = nanmean(diff,1)' .* h.data;
h.data = diff;
afni_niml_writesimple(h,['mask_p05_group_mean_diff_' out]);

out = out_rh;
h = afni_niml_readsimple(out_rh);
h.data(abs(h.data) <= cutoff) = 0;
h.data(h.data ~= 0) = 1;
afni_niml_writesimple(h,['mask_p05_' out]);
diff = surf_ds_rh.samples(1:29,:) - surf_ds_rh.samples(30:58,:);
diff = nanmean(diff,1)' .* h.data;
h.data = diff;
afni_niml_writesimple(h,['mask_p05_group_mean_diff_' out]);




