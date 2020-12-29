%% Get the homotopic node-to-node correspondences between hemispheres
% Inspired by method from Jo et al. (2012) https://doi.org/10.1371/journal.pone.0048847 
% Left-Right node correspondence based on most similar distance profile to
% all FS ROI centers. 

% Note: This takes a long time and scales with mesh density. Over 8.5 hours
% for just ONE hemipshere of the ld141 mesh on our iMac with 8-core 3.6 GHz
% Intel Core i9, 32GB.


purge

fs_dir = '/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs';
cd(fs_dir)

hemi = {'lh' 'rh'};
dens = {'141'}; % '60'

out_dir = [fs_dir '/FS_node_to_ROI_distances'];
% mkdir(out_dir);

for dd = 1:numel(dens)
for hh = 1:numel(hemi)
    d = dens{dd};
    h = hemi{hh};
    C = readtable(['FS_ROI_centers_' h '_' d '.txt']);
    prefix = [out_dir '/' h '_ld' d];
    
    % Load example surface dset
    if strcmp(d,'60')
        dset = afni_niml_readsimple(['GroupMean_stats_' h '_REML_Coeff_Dp.niml.dset']);
    elseif strcmp(d,'141')
        dset = afni_niml_readsimple(['GroupMean_tracks_ss3t_50M.wholebrain_length_map.al2anat.' h '.6mm.niml.dset']);
    end
    nodes = dset.node_indices;
    
    % Loop through each ROI, calculate distance from every surface node
    parfor nn = 1:height(C)
        % Specify the node pairs
        center = C.Var4(nn);
        roi_ind = num2str(nn,'%02.f');
        roi_label = C.Var1{nn};
        pairs = [repelem(nodes,numel(center)), repmat(center,numel(nodes),1)];
        fname = [prefix '_ROI_' roi_ind '_' roi_label];
        writematrix(pairs,[fname '_node_pairs.txt']);
        unix(['SurfDist -i std.' d '.' h '.smoothwm.gii -input ' fname '_node_pairs.txt > ' fname '_node_dists.txt']);
    end
end
end


%% Now load the distances and find distance similarities
dens = {'ld60'};%,'ld141'}
cd(out_dir)
for dd = 1:numel(dens)
    d = dens{dd};    
    % Load surface dsets to use as template
    if strcmp(d,'ld60')
        dset_lh = afni_niml_readsimple([fs_dir '/GroupMean_stats_lh_REML_Coeff_Dp.niml.dset']);
        dset_rh = afni_niml_readsimple([fs_dir '/GroupMean_stats_rh_REML_Coeff_Dp.niml.dset']);
    elseif strcmp(d,'ld141')
        dset_lh = afni_niml_readsimple([fs_dir '/GroupMean_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset']);
        dset_rh = afni_niml_readsimple([fs_dir '/GroupMean_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset']);
    end
    
    % Load the distance data from the individual ROI files
    lh_dists = dir(['lh_' d '*node_dists.txt']);
    rh_dists = dir(['rh_' d '*node_dists.txt']);
    for ii = 1:numel(lh_dists)
        l = readmatrix(lh_dists(ii).name);
        r = readmatrix(rh_dists(ii).name);
        l_vecs(:,ii) = l(:,3);
        r_vecs(:,ii) = r(:,3);
    end
    
    % Run Pearson correlation across vectors, creates an NxN matrix (N =
    % # of nodes)
    lr_similarity = corr(l_vecs',r_vecs');
    
    % Find max rho in the left hemisphere for each node in right
    for ii = 1:size(lr_similarity,1)
        node_pairs_LtoR(ii,1) = ii;
        node_pairs_RtoL(ii,2) = ii;
        [~,node_pairs_LtoR(ii,2)] = max(lr_similarity(ii,:));
        [~,node_pairs_RtoL(ii,1)] = max(lr_similarity(:,ii));
    end
    
    % Populate the datasets
    dset_lh.data = (1:length(dset_lh.node_indices))';
    dset_rh.data = node_pairs_LtoR(:,2);
    % Write dsets based on homotopic correspondence. 
    afni_niml_writesimple(dset_lh,['homotopic_correspondence_LtoR_' d '_lh.niml.dset']);
    afni_niml_writesimple(dset_rh,['homotopic_correspondence_LtoR_' d '_rh.niml.dset']);

    % Populate the datasets
    dset_rh.data = (1:length(dset_rh.node_indices))';
    dset_lh.data = node_pairs_RtoL(:,1);
    % Write dsets based on homotopic correspondence. 
    afni_niml_writesimple(dset_lh,['homotopic_correspondence_RtoL_' d '_lh.niml.dset']);
    afni_niml_writesimple(dset_rh,['homotopic_correspondence_RtoL_' d '_rh.niml.dset']);

    % Save node mappings (subtract 1 so the values correspond to 0-based
    % node index used by AFNI/SUMA)
    writematrix(node_pairs_LtoR-1,['homotopic_correspondence_LtoR_' d '.txt']);
    writematrix(node_pairs_RtoL-1,['homotopic_correspondence_RtoL_' d '.txt']);
end
    


%% Map data from one surface to the other
% NOTE - The order of the process depends

clear all
cd('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/AllSubs_dsets')

% Load mapping files and add back 1 to account for 0-based index
Lseed_Rtarg = readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_LtoR_ld60.txt');
Ltarg_Rseed = readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_RtoL_ld60.txt');
Lseed_Rtarg = Lseed_Rtarg+1;
Ltarg_Rseed = Ltarg_Rseed+1;



% ----- RH data to LH surface -----
data = {'AllSubs_std.60.rh.PP19_Dp-Da_math.MNI152.votc.inflated.14mm_diam.niml.dset'
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da.niml.dset'
        'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp.niml.dset'
        'homotopic_correspondence_RtoL_ld60_rh.niml.dset'};
for ii = 1:numel(data)
    d = afni_niml_readsimple(data{ii});
    dnew = d;
    for kk = 1:size(d.data,2)
    for jj = 1:numel(d.node_indices)
       % Find all the seed nodes that mapped to this target
       inds = find(Ltarg_Rseed(:,1) == jj);
       if numel(inds) > 0
           val = mean(d.data(inds,kk));
       else
       % If there were no seeds mapped to this target, base the data on
       % this node's target
           val = d.data(Lseed_Rtarg(jj,2),kk);
       end
       dnew.data(jj,kk) = val;
    end
    end
    afni_niml_writesimple(dnew,strrep(data{ii},'.niml.dset','_MAPPED2CONTRA.niml.dset'));
end



% ----- LH data to RH surface -----
data = {'AllSubs_std.60.lh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.niml.dset'
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset'
        'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset'
        'homotopic_correspondence_LtoR_ld60_lh.niml.dset'};
for ii = 1:numel(data)
    d = afni_niml_readsimple(data{ii});
    dnew = d;
    for kk = 1:size(d.data,2)
    for jj = 1:numel(d.node_indices)
       % Find all the seed nodes that mapped to this target
       inds = find(Lseed_Rtarg(:,2) == jj);
       if numel(inds) > 0
           val = mean(d.data(inds,kk));
       else
       % If there were no seeds mapped to this target, base the data on
       % this node's target
           val = d.data(Ltarg_Rseed(jj,1),kk);
       end
       dnew.data(jj,kk) = val;
    end
    end
    afni_niml_writesimple(dnew,strrep(data{ii},'.niml.dset','_MAPPED2CONTRA.niml.dset'));
end





%% Verification
% Confirmed that compared to previous iterations, the final method (v3)
% yielded transformed distributions that were closest to original data
%
% v0 = afni_niml_readsimple('AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp.niml.dset');
% v1 = afni_niml_readsimple('AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp_MAPPED2CONTRA.niml.dset');
% v2 = afni_niml_readsimple('AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp_MAPPED2CONTRA_v2.niml.dset');
% v3 = afni_niml_readsimple('AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp_MAPPED2CONTRA_v3.niml.dset');
% 
% for ii = 1:29
%     [~,p(ii,1),kstat(ii,1)] = kstest2(v0.data(:,ii),v1.data(:,ii));
%     [~,p(ii,2),kstat(ii,2)] = kstest2(v0.data(:,ii),v2.data(:,ii));
%     [~,p(ii,3),kstat(ii,3)] = kstest2(v0.data(:,ii),v3.data(:,ii));
% end























    
    
    
    