%% Make circular ROIs at calculated FS region centers for verification

out_dir = '/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs';
cd(out_dir);
% LH
T = readtable('FS_ROI_centers_lh_60.txt');
writematrix(T.Var4,'FS_lh_center_nodes.txt');
writematrix(ones(height(T),1),'FS_lh_center_node_tmpval.txt');
% RH
T = readtable('FS_ROI_centers_rh_60.txt');
writematrix(T.Var4,'FS_rh_center_nodes.txt');
writematrix(ones(height(T),1),'FS_rh_center_node_tmpval.txt');


hemi = 'lh';
out = 'FS_lh_centers';
out_dir = '/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs';
MNI_dir = '/Users/benconrad/.afni/data/suma_MNI152_2009';
mesh = '60';
sid = 'MNI152_2009';
spec = ['std.' mesh '.' sid '_' hemi '.spec']; % Surface spec file
surf = ['std.' mesh '.' hemi '.inflated.gii']; % Surface to "grow" ROIs on
d = '5';
cd(MNI_dir);
unix(['ROIgrow -overwrite'...
    ' -insphere ' d...
    ' -spec ' spec ...
    ' -surf ' surf ...
    ' -roi_nodes ' out_dir '/FS_lh_center_nodes.txt'...
    ' -roi_labels ' out_dir '/FS_lh_center_node_tmpval.txt'...
    ' -prefix ' out_dir '/' out '.inflated.' d 'mm_diam']);


hemi = 'rh';
out = 'FS_rh_centers';
spec = ['std.' mesh '.' sid '_' hemi '.spec']; % Surface spec file
surf = ['std.' mesh '.' hemi '.inflated.gii']; % Surface to "grow" ROIs on
unix(['ROIgrow -overwrite'...
    ' -insphere ' d...
    ' -spec ' spec ...
    ' -surf ' surf ...
    ' -roi_nodes ' out_dir '/FS_rh_center_nodes.txt'...
    ' -roi_labels ' out_dir '/FS_rh_center_node_tmpval.txt'...
    ' -prefix ' out_dir '/' out '.inflated.' d 'mm_diam']);


