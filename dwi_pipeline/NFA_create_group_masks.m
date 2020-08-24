%% Group mask creation 
% Based on greater than 50 percent of subjects having at least one streamline to surface node
% Then removal of the seed ROIs 

roi_strings = {'1*/tracks_ss3t_*_combined_Dp-Da.lh.TDI_ends.norm.al2anat.lh.niml.dset'
               '1*/tracks_ss3t_*_combined_Dp-Da.lh.TDI_ends.norm.al2anat.rh.niml.dset'
               '1*/tracks_ss3t_*_combined_Lp-La.lh.TDI_ends.norm.al2anat.lh.niml.dset'
               '1*/tracks_ss3t_*_combined_Lp-La.lh.TDI_ends.norm.al2anat.rh.niml.dset'
               '1*/tracks_ss3t_*_combined_Dp-Da.rh.TDI_ends.norm.al2anat.lh.niml.dset'
               '1*/tracks_ss3t_*_combined_Dp-Da.rh.TDI_ends.norm.al2anat.rh.niml.dset'
               '1*/tracks_ss3t_*_combined_Lp-La.rh.TDI_ends.norm.al2anat.lh.niml.dset'
               '1*/tracks_ss3t_*_combined_Lp-La.rh.TDI_ends.norm.al2anat.rh.niml.dset'};
           
labels = {'Dp-Da.lh.lh.niml.dset'
           'Dp-Da.lh.rh.niml.dset'
           'Lp-La.lh.lh.niml.dset'
           'Lp-La.lh.rh.niml.dset'
           'Dp-Da.rh.lh.niml.dset'
           'Dp-Da.rh.rh.niml.dset'
           'Lp-La.rh.lh.niml.dset'
           'Lp-La.rh.rh.niml.dset'};

% Get non-zero count and greater than 50% nonzero masks for each seed
% ROI/hemisphere (n = 29 subjects so 15 is > 50%)
for ii = 1:numel(roi_strings)
    unix(['3dMean -count -prefix group_nzero_count.' labels{ii} ' ' roi_strings{ii}])
    unix(['3dcalc -prefix group_50perc_nzero.' labels{ii} ' -a group_nzero_count.' labels{ii} ' -expr "ispositive(a-14)"'])
end

unix(['3dmean -mask_union -prefix group_50perc_union.lh.lh.niml.dset'...
     ' group_50perc_nzero.Dp-Da.lh.lh.niml.dset'...
     ' group_50perc_nzero.Lp-La.lh.lh.niml.dset']);
 
unix(['3dmean -mask_union -prefix group_50perc_union.lh.rh.niml.dset'...
     ' group_50perc_nzero.Dp-Da.lh.rh.niml.dset'...
     ' group_50perc_nzero.Lp-La.lh.rh.niml.dset']);
 
 unix(['3dmean -mask_union -prefix group_50perc_union.rh.rh.niml.dset'...
     ' group_50perc_nzero.Dp-Da.rh.rh.niml.dset'...
     ' group_50perc_nzero.Lp-La.rh.rh.niml.dset']);
 
 unix(['3dmean -mask_union -prefix group_50perc_union.rh.lh.niml.dset'...
     ' group_50perc_nzero.Dp-Da.rh.lh.niml.dset'...
     ' group_50perc_nzero.Lp-La.rh.lh.niml.dset']);
 
 
 
 unix(['3dmean -mask_union -prefix group_50perc_union.lh.niml.dset'...
     ' group_50perc_nzero.Dp-Da.lh.lh.niml.dset'...
     ' group_50perc_nzero.Lp-La.lh.lh.niml.dset'...
     ' group_50perc_nzero.Dp-Da.rh.lh.niml.dset'...
     ' group_50perc_nzero.Lp-La.rh.lh.niml.dset']);
 
  
 unix(['3dmean -mask_union -prefix group_50perc_union.rh.niml.dset'...
     ' group_50perc_nzero.Dp-Da.lh.rh.niml.dset'...
     ' group_50perc_nzero.Lp-La.lh.rh.niml.dset'...
     ' group_50perc_nzero.Dp-Da.rh.rh.niml.dset'...
     ' group_50perc_nzero.Lp-La.rh.rh.niml.dset']);
 
 
%% Create union of subject ROI masks and remove from 50% maps created above
% LH
roi_files = dir('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/*_proc/*.freesurfer/SUMA/std.141.lh.PollackandPrice19_*p-*a.MNI152.votc.inflated.14mm_diam.1.1D');
nodes = [];
for ii = 1:numel(roi_files)
    fid = fopen([roi_files(ii).folder '/' roi_files(ii).name]);
    cdata = textscan(fid,'%f%f','delimiter',',', 'HeaderLines',2);
    fclose(fid);
    % Get indices for this subject
    nodes = [nodes;cdata{1,1}]; % Add 1, as these indices start at 0, whereas stats starts at 1...
end   
nodes_unique = unique(nodes);
dlmwrite('lh_rois.1D',[nodes_unique,ones(size(nodes_unique))],'Delimiter','\t','precision',10)
unix(['ConvertDset -o_niml -pad_to_node ld141'...
      ' -input lh_rois.1D -prefix lh_rois -node_index_1D lh_rois.1D[0]']);

unix('3dcalc -prefix lh_final_group_mask.niml.dset -a group_50perc_union.lh.niml.dset -b lh_rois.niml.dset -expr "a*not(b)"');


% RH
roi_files = dir('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/*_proc/*.freesurfer/SUMA/std.141.rh.PollackandPrice19_*p-*a.MNI152.votc.inflated.14mm_diam.1.1D');
nodes = [];
for ii = 1:numel(roi_files)
    fid = fopen([roi_files(ii).folder '/' roi_files(ii).name]);
    cdata = textscan(fid,'%f%f','delimiter',',', 'HeaderLines',2);
    fclose(fid);
    % Get indices for this subject
    nodes = [nodes;cdata{1,1}]; % Add 1, as these indices start at 0, whereas stats starts at 1...
end   
nodes_unique = unique(nodes);
dlmwrite('rh_rois.1D',[nodes_unique,ones(size(nodes_unique))],'Delimiter','\t','precision',10)
unix(['ConvertDset -o_niml -pad_to_node ld141'...
      ' -input rh_rois.1D -prefix rh_rois -node_index_1D rh_rois.1D[0]']);

unix('3dcalc -prefix rh_final_group_mask.niml.dset -a group_50perc_union.rh.niml.dset -b rh_rois.niml.dset -expr "a*not(b)"');

