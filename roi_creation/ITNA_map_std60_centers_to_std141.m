%% Map std60 mesh centers to std141 mesh
% To avoid having to run the distance analysis in std141 space (would take
% a very long time)

purge
hemi = {'lh' 'rh'};
cd('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs')

for hh = 2%1:numel(hemi)
    h = hemi{hh};
    % Create sorted list of nodes/ROIs from std60 mesh
    C = readtable(['FS_ROI_centers_' h '_60.txt']);
    C = C.Var4;
    C = [C,(1:length(C))'];
    [~,I] = sort(C(:,1));
    C = C(I,:);
    writematrix(C(:,1),['FS_' h '_center_nodes_sort.txt']);
    writematrix(C,['FS_' h '_center_nodes_sort_rois.txt']);
    
    % Convert to full node list dataset
    unix(['ConvertDset -o_1D -pad_to_node ld60 '...
        ' -input FS_' h '_center_nodes_sort_rois.txt '...
        ' -node_index_1D FS_' h '_center_nodes_sort.txt '...
        ' -prefix FS_' h '_center_nodes_allrows']);
    
    % Perform mapping to surface
    unix(['SurfToSurf -i std.141.' h '.inflated.gii -sv MNI152_2009_SurfVol.nii '...
        '             -i std.60.' h '.inflated.gii -sv MNI152_2009_SurfVol.nii'...
        ' -data  FS_' h '_center_nodes_allrows.1D.dset'...
        ' -prefix mapping_std60_to_std141_' h...
        ' -output_params NearestNode']);
    
    % Get the node index corresponding to each ROI
    T = readmatrix(['mapping_std60_to_std141_' h '.1D'],'FileType','text','NumHeaderLines',17);
    roi_list = 1:74;
    for ii = 1:numel(roi_list)
        out(ii,1) = find(T(:,4) == ii, 1, 'first');
        out(ii,2) = ii;
    end
    
    % Write out info to new text files
    writematrix(out(:,1),['FS_' h '_center_nodes_141.txt']);
    C = readtable(['FS_ROI_centers_' h '_60.txt']);
    C.Var4 = out(:,1);
    writetable(C,['FS_ROI_centers_' h '_141.txt']);
end




            