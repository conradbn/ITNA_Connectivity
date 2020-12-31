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
dens = {'ld141'}; %'60'
cd(out_dir)
for dd = 1:numel(dens)
    d = dens{dd};    
    % Load example surface dsets to use as template
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
    
    % Make single precision to reduce memory load
    l_vecs = single(l_vecs);
    r_vecs = single(r_vecs);
    
    if strcmp(d,'ld60')
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
    % For high density mesh, need to go one node at a time (i.e., can't hold a
    % 199k x 199k correlation matrix in memory)
    elseif strcmp(d,'ld141')
        parfor ii = 1:size(l_vecs,1)
            disp(['Working on node # ' num2str(ii)]);
            col1_LtoR(ii,1) = ii;
            col2_RtoL(ii,1) = ii;
            nodeL = l_vecs(ii,:);
            nodeR = r_vecs(ii,:);
            corrLnode = corr(nodeL',r_vecs');
            corrRnode = corr(nodeR',l_vecs');
            [~,ind] = max(corrLnode);
            [~,col2_LtoR(ii,1)] = max(corrLnode);
            [~,col1_RtoL(ii,1)] = max(corrRnode);
        end
        node_pairs_LtoR = [col1_LtoR,col2_LtoR];
        node_pairs_RtoL = [col1_RtoL,col2_RtoL];
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

%% Map data to contralateral hemishperes
clear all
cd('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/AllSubs_dsets')

% Set data files to map, including mesh density and mapping direction
data = {'ld60'  'RtoL' 'AllSubs_std.60.rh.PP19_Dp-Da_math.MNI152.votc.inflated.14mm_diam.niml.dset'
        'ld60'  'RtoL' 'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zdiff.Dp-Da.niml.dset'
        'ld60'  'RtoL' 'AllSubs_Dp-Da_math.rh.beta_series_corr.rh.Zmap.Dp.niml.dset'
        'ld60'  'RtoL' 'homotopic_correspondence_RtoL_ld60_rh.niml.dset'
        'ld60'  'LtoR' 'AllSubs_std.60.lh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.niml.dset'
        'ld60'  'LtoR' 'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zdiff.Dp-Da.niml.dset'
        'ld60'  'LtoR' 'AllSubs_Dp-Da.lh.beta_series_corr.lh.Zmap.Dp.niml.dset'
        'ld60'  'LtoR' 'homotopic_correspondence_LtoR_ld60_lh.niml.dset'
        'ld141' 'RtoL' 'Allsubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
        'ld141' 'RtoL' 'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c.niml.dset'
        'ld141' 'RtoL' 'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset'
        'ld141' 'RtoL' 'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.niml.dset'
        'ld141' 'RtoL' 'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.log.niml.dset'
        'ld141' 'LtoR' 'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
        'ld141' 'LtoR' 'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.log+c.niml.dset'
        'ld141' 'LtoR' 'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset'
        'ld141' 'LtoR' 'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.niml.dset'
        'ld141' 'LtoR' 'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.log.niml.dset'};
    
% Loop through each dataset
for ii = 1:numel(data)
    % Load data structure
    data_struct = afni_niml_readsimple(data{ii,3});
    % Get mapping direction
    direction = data{ii,2};
    % Check which mesh density and load mapping files (add back 1 to
    % account for 0-based index)
    if strcmp(data{ii,1},'ld60')
        Lseed_Rtarg = 1 + readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_LtoR_ld60.txt');
        Ltarg_Rseed = 1 + readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_RtoL_ld60.txt');
    elseif strcmp(data{ii,1},'ld141')
        Lseed_Rtarg = 1 + readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_LtoR_ld141.txt');
        Ltarg_Rseed = 1 + readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_RtoL_ld141.txt');
    end
    % Run the mapping function
    data_struct.data = single(map_data(data_struct,Ltarg_Rseed,Lseed_Rtarg,direction));
    % Write new data file
    afni_niml_writesimple(data_struct,strrep(data{ii,3},'.niml.dset','_MAP2CON.niml.dset'));
end

%% Mapping function
function mapped_data = map_data(data_struct,LtRs,LsRt,direction)
d = data_struct.data;
mapped_data = zeros(size(d,1),size(d,2));
for ss = 1:size(d,2) % Subjects/volumes
    disp(['Mapping data (direction ' direction ') for Subject/Volume # ' num2str(ss)])
    ds = d(:,ss);
    mapped = zeros(size(d,1),1);
    parfor nn = 1:size(d,1) % Nodes
        if strcmp(direction,'RtoL')
            % Find all the seed nodes that mapped to this target
            inds = find(LtRs(:,1) == nn);
            if numel(inds) > 0
                val = mean(ds(inds));
            else
                % If there were no seeds mapped to this target, base the data on
                % this node's target
                val = ds(LsRt(nn,2));
            end
        elseif strcmp(direction,'LtoR')
            % Find all the seed nodes that mapped to this target
            inds = find(LsRt(:,2) == nn);
            if numel(inds) > 0
                val = mean(ds(inds));
            else
                % If there were no seeds mapped to this target, base the data on
                % this node's target
                val = ds(LtRs(nn,1));
            end
        end
        mapped(nn,1) = val;
    end
    mapped_data(:,ss) = mapped;
end
end
              