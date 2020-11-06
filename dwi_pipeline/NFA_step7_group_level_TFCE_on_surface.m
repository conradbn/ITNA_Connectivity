purge
hemi = {'lh','rh'};

%% CONTRASTS
for ii = 1:numel(hemi)
    h = hemi{ii};
    % Create combined, covariate corrected residual dataset
    cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/TFCE')
    out = ['allsubs_' h '_Digit_Letter.' h '.TDI_ends.norm.6mm.log+c.niml.dset'];
    unix(['3dTcat -overwrite -prefix ' out ...
         ' ../*/*/resid_from_cov_Digit_TDI_ends_surf_' h '_6mm_log+c_3dttest++.' h '.mask.niml.dset'...
         ' ../*/*/resid_from_cov_Letter_TDI_ends_surf_' h '_6mm_log+c_3dttest++.' h '.mask.niml.dset']);

    % Load in the subject data
    surf_ds = afni_niml_readsimple(out);
    surf_ds.labels = num2cell(1:size(surf_ds.data,2));
    surf_ds = cosmo_surface_dataset(surf_ds);

    % Get faces and vertices info from MNI surface gii file
    [vertices,faces]=surfing_read(['../suma_MNI152_2009/std.141.' h '.smoothwm.gii']);
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
    out = ['tfce_' h '_D-L_log+c.' h '.1k.noroimasks.niml.dset'];
    niters = 1000;

    % Mask the dataset to reduce to only areas of "consistent" connectivity to
    % digit OR letter ROIs
    mask_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Data';
    m1 = [mask_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_Dp-Da.' h '.TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
    m2 = [mask_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_Lp-La.' h '.TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
    m3 = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Masks/' h '_rois.niml.dset'];

    m1s = afni_niml_readsimple(m1);
    m2s = afni_niml_readsimple(m2);
    m3s = afni_niml_readsimple(m3);
    ind_keep = m1s.data == 1 | m2s.data == 1;
%     ind_keep(m3s.data(:,2)==1) = 0;
    surf_ds.samples(:,~ind_keep) = NaN;

    % Run TFCE
    NFA_run_TFCE(surf_ds,vertices,faces,niters,out,'');
    
    % Create masks of significant difference
    cutoff = 2.5758;%1.96;
    h = afni_niml_readsimple(out);
    h.data(abs(h.data) <= cutoff) = 0;
    h.data(h.data ~= 0) = 1;
    afni_niml_writesimple(h,['mask_p01_' out]);
    diff = surf_ds.samples(1:29,:) - surf_ds.samples(30:58,:);
    diff = nanmean(diff,1)' .* h.data;
    h.data = diff;
    afni_niml_writesimple(h,['mask_p01_group_mean_diff_' out]);
end

%% CONJUNCTIONS
% for ii = 1:numel(hemi)
%     h = hemi{ii};
%     % Create combined, covariate corrected residual dataset
%     cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/TFCE')
%     out = ['allsubs_' h '_Digit_Letter_CONJ.' h '.TDI_ends.norm.6mm.log+c.niml.dset'];
%     unix(['3dTcat -overwrite -prefix ' out ...
%          ' ../*/tracks_ss3t_*_50M_Dp-Da.' h '.TDI_ends.norm.al2anat.' h '.6mm.log+c.niml.dset'...
%          ' ../*/tracks_ss3t_*_50M_Lp-La.' h '.TDI_ends.norm.al2anat.' h '.6mm.log+c.niml.dset']);
%     
%     
%     % Load in the subject data, and create thresholded conjunction map for
%     % each subject (binary)
%     t = -7;
%     surf_ds = afni_niml_readsimple(out);
%     for jj = 1:29
%         s1 = surf_ds.data(:,jj);
%         s2 = surf_ds.data(:,jj*2);
%         surf_ds_bin(:,jj) = s1 > t & s2 > t;
%     end
%     
%     surf_ds.data = surf_ds_bin;
%     out_bin = strrep(out,'CONJ','CONJ_BIN');
%     afni_niml_writesimple(surf_ds,out_bin);
%     
%     surf_ds.labels = num2cell(1:size(surf_ds.data,2));
%     surf_ds = cosmo_surface_dataset(surf_ds);
% 
%     % Get faces and vertices info from MNI surface gii file
%     [vertices,faces]=surfing_read(['../suma_MNI152_2009/std.141.' h '.smoothwm.gii']);
%     faces = double(faces);
%     vertices = double(vertices);
% 
%     % Set up surface dataset structure, to prepare it for Cosmo TFCE
%     nsubs = size(surf_ds.samples,1);
%     surf_ds.fa.center_ids = surf_ds.fa.node_indices;
%     surf_ds.sa.chunks = (1:nsubs)';
%     surf_ds.sa.targets = ones(nsubs,1);
% 
%     % Set filename for the output statistic dataset and number of iterations
%     out = ['tfce_' h '_D+L_log+c.' h '.100k.niml.dset'];
%     niters = 100000;
% 
%     % Mask the dataset to reduce to only areas of "consistent" connectivity to
%     % digit OR letter ROIs
%     mask_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Data';
%     m1 = [mask_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_Dp-Da.' h '.TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
%     m2 = [mask_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_Lp-La.' h '.TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
%     m3 = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Masks/' h '_rois.niml.dset'];
% 
%     m1s = afni_niml_readsimple(m1);
%     m2s = afni_niml_readsimple(m2);
%     m3s = afni_niml_readsimple(m3);
%     ind_keep = m1s.data == 1 & m2s.data == 1; % Where both are true
%     ind_keep(m3s.data(:,2)==1) = 0;
%     surf_ds.samples(:,~ind_keep) = NaN;
%     
%     % Run TFCE
%     NFA_run_TFCE(surf_ds,vertices,faces,niters,out,'');
%     
%     % Create masks of significant Conjunction
%     cutoff = norminv(0.995);%1.96;
%     h = afni_niml_readsimple(out);
%     h.data(abs(h.data) <= cutoff) = 0;
%     h.data(h.data ~= 0) = 1;
%     afni_niml_writesimple(h,['mask_p01_' out]);
% end


