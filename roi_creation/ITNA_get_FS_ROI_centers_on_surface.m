%% Get FreeSurfer region centers on MNI152-2009 surfaces
% This script finds the central node of every FS region (aparc.a2009s parcellation)
% for both hemispheres and for both standard mesh densities (ld60 and
% ld141). The purpose is to use these central nodes as landmarks for
% determining the homotopic node correspondance between hemispheres, as in
% Jo et al. (2012) https://doi.org/10.1371/journal.pone.0048847 

purge;

% Get unique region names
lut = readmatrix('aparc.a2009s+aseg_REN_all.niml.lt','FileType','text','Range','6:202','OutputType','char');
ctx = lut(contains(lut(:,2),'ctx_'),2);
ctx_rois = cellfun(@(S) S(8:end), ctx, 'Uniform', 0);
ctx_rois = unique(ctx_rois);
ctx_rois(contains(ctx_rois,'Unknown')) = [];

hemi = {'lh','rh'}; %{'lh','rh'};
dens = {'60'}; %{'60','141'};


for hh = 1:numel(hemi)
for dd = 1:numel(dens)
    h = hemi{hh};
    d = dens{dd};
    
    % Get SUMA MNI FreeSurfer Parcellation information
    FS = afni_niml_read(['std.' d '.' h '.aparc.a2009s.annot.niml.dset']);
    FS_rois = FS.nodes{1}.data;
    FS_nodes = FS.nodes{2}.data;
    ctx_rois_ids = readmatrix('FS_idcodes.csv','OutputType','char');
    
    % Build table of all the ROI info/center node
    center_node_table = table();
    parfor ii = 1:numel(ctx_rois)
        
        disp(['hemi-' h ' density-' d ' roi # ' num2str(ii) ' of ' num2str(numel(ctx_rois))]);
        % Load the node information
        c = ctx_rois{ii};
        roi_label = ['ctx_' h '_' c];
        id_renum = str2double(lut{strcmp(lut(:,2),roi_label),1});
        id = ctx_rois_ids(contains(ctx_rois_ids(:,2),[h '_' c]),1);
        id = str2double(id);

        % Get the nodes from that ROI and find all uniqure pairings
        roi_nodes = FS_nodes(FS_rois == id);
        out=nchoosek(roi_nodes,2);
        writematrix(out,['tmp_node_pairs_' num2str(ii) '.txt']);

        % Get distance for every pair
        unix(['SurfDist -i std.' d '.' h '.smoothwm.gii -input tmp_node_pairs_' num2str(ii) '.txt > tmp_dists_' num2str(ii) '.1D']);

        % Load distances
        dists = readmatrix(['tmp_dists_' num2str(ii) '.1D'],'FileType','text');
        dist_sum = [];
        for jj = 1:numel(roi_nodes)
            r = roi_nodes(jj);
            rind = any((dists(:,1:2) == r),2);
            dist_sum(jj) = sum(dists(rind,3));
        end

        % Find the node with the smallest total distance
        [~,min_dist] = min(dist_sum);
        center_node = roi_nodes(min_dist);

        % Populate info into temporary table
        tmpT = table();
        tmpT.label = {roi_label};
        tmpT.id = id;
        tmpT.id_renum = id_renum;
        tmpT.center = center_node;

        % Add to full table
        center_node_table(ii,:) = tmpT;
    end
    writetable(center_node_table,['FS_ROI_centers_' h '_' d '.txt'],'Delimiter','tab');
    unix('rm -f tmp*');
end
end