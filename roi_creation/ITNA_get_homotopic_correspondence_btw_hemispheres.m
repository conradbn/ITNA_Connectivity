%% Get the homotopic node-to-node correspondences between hemispheres
% Inspired by method from Jo et al. (2012) https://doi.org/10.1371/journal.pone.0048847 
% Left-Right node correspondence based on most similar distance profile to
% all FS ROI centers. 


purge

fs_dir = '/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs';
cd(fs_dir)

fs_centers_lh = readtable('FS_ROI_centers_lh_60.txt'); 
fs_centers_rh = readtable('FS_ROI_centers_rh_60.txt'); 

hemi = {'lh' 'rh'};
dens = {'60'};

for dd = 1:numel(dens)
for hh = 1:numel(hemi)
    d = dens{dd};
    h = hemi{hh};
    C = readtable(['FS_ROI_centers_' h '_60.txt']);
    dset = afni_niml_readsimple(['GroupMean_stats_' h '_REML_Coeff_Dp.niml.dset']); % Load example surface dset
    % Specify the node pairs
    nodes = dset.node_indices;
    centers = C.Var4;
    pairs = [repelem(nodes,numel(centers)), repmat(centers,numel(nodes),1)];
    writematrix(pairs,'tmp_node_pairs.txt');
    unix(['SurfDist -i std.' d '.' h '.smoothwm.gii -input tmp_node_pairs.txt > FS_ROI_center_dists_' d '_' h '.txt']);
end
end


