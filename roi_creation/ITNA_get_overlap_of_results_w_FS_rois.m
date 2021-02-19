purge;

% Get the names of all TFCE files
tfce_maps = dir('/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Results/*/*TFCE_Zscore_100000iters_MASK.niml.dset');

% Get unique region names from SUMA MNI FreeSurfer Parcellation
cd('/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/FreeSurfer_ROIs')
lut = readmatrix('aparc.a2009s+aseg_REN_all.niml.lt','FileType','text','Range','6:202','OutputType','char');
ctx = lut(contains(lut(:,2),'ctx_'),2);
ctx_rois = cellfun(@(S) S(8:end), ctx, 'Uniform', 0);
ctx_rois = unique(ctx_rois);
ctx_rois(contains(ctx_rois,'Unknown')) = [];
%%
for tt =34% 1:numel(tfce_maps)
    zscr_map = afni_niml_readsimple([tfce_maps(tt).folder '/' tfce_maps(tt).name]);
    % Set density and hemisphere
    if size(zscr_map.data,1) == 36002
        d = '60';
    else
        d = '141';
    end
    if contains(tfce_maps(tt).name,'on_RHsurf')
        h = 'rh';
    else
        h = 'lh';
    end
   
    
    % Get SUMA MNI FreeSurfer Parcellation information
    FS = afni_niml_read(['std.' d '.' h '.aparc.a2009s.annot.niml.dset']);
    FS_rois = FS.nodes{1}.data;
    FS_nodes = FS.nodes{2}.data;
    ctx_rois_ids = readmatrix('FS_idcodes.csv','OutputType','char');
    
    % For each ROI, load the data
    T = table(  );
    for ii = 1:numel(ctx_rois)
        
        disp(['hemi-' h ' density-' d ' roi # ' num2str(ii) ' of ' num2str(numel(ctx_rois))]);
        % Load the node information
        c = ctx_rois{ii};
        roi_label = ['ctx_' h '_' c];
        id_renum = str2double(lut{strcmp(lut(:,2),roi_label),1});
        id = ctx_rois_ids(contains(ctx_rois_ids(:,2),[h '_' c]),1);
        id = str2double(id);

        % Get the nodes from that ROI and find all uniqure pairings
        roi_nodes = FS_rois == id;
        
        % Add stats to the table
        T.region(ii) = string(roi_label);
        T.count_all(ii) = sum(roi_nodes);
        T.count_pos(ii) = sum(zscr_map.data(roi_nodes) > 1.96);
        T.count_neg(ii) = sum(zscr_map.data(roi_nodes) < -1.96);
        T.perc_pos(ii) = 100 * T.count_pos(ii)/T.count_all(ii);
        T.perc_neg(ii) = 100 * T.count_neg(ii)/T.count_all(ii);
    end
end