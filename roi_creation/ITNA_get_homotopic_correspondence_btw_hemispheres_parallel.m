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

out_dir = 'FS_node_to_ROI_distances';
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
dens = {'ld60','ld141'};
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
    
    % Find max rho for each node
    for ii = 1:size(lr_similarity,1)
        lr_node_pairs(ii,1) = ii;
        [~,lr_node_pairs(ii,2)] = max(lr_similarity(:,ii));
    end

    % Populate RH dataset based on corresponding LH node values
    dset_lh.data = (1:length(dset_lh.node_indices))';
    dset_rh.data(:) = 0;
    for ii = 1:length(dset_lh.node_indices)
        l = lr_node_pairs(ii,1);
        r = lr_node_pairs(ii,2);
        dset_rh.data(r,1) = l;
    end
    
    % Write dsets based on homotopic correspondence. 
    afni_niml_writesimple(dset_lh,['homotopic_correspondence_' d '_lh_.niml.dset']);
    afni_niml_writesimple(dset_rh,['homotopic_correspondence_' d '_rh_.niml.dset']);
end
    
    
    
    
    
    
    