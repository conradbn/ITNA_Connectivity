%% Use ROIgrow command to create circular ROIs on the cortical surfaces from FreeSurfer
% Center on standard surface (std.60 density) nodes based on Pollack &
% Price 2019, NeuroImage
% This script uses the courser grained "std.60" mesh density, which is that
% used for fMRI analyses. Ultimately need these masks for the gPPI analysis.

purge

DpDa_math_60 = '/Users/nbl_imac2/Documents/GitHub/NFA_Stucture_Function/roi_creation/rois/LitCoord_Digit_Pollack19_54_-52_-14_std.60_rh';
DpDa_math_141 = '/Users/nbl_imac2/Documents/GitHub/NFA_Stucture_Function/roi_creation/rois/LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh';
diam = '14';

cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData');
sub_dirs = dir('*_proc');

% Loop through subjects
for ii = 1:numel(sub_dirs)
    % Get subject ID and go to freesurfer SUMA folder
    sid = strsplit(sub_dirs(ii).name,'_'); sid = sid{1,1};
    disp(['** Working on subject ' sid '... **']);
    cd(['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sid '_proc/' sid '.freesurfer/SUMA']); 
    
    % Create circular ROI on the surface, around selected node in each
    % hemisphere
    spec = ['std.60.' sid '_rh.spec']; % Surface spec file
    surf = 'std.60.rh.inflated.gii'; % Surface to "grow" ROIs on
    unix(['ROIgrow -overwrite'... 
          ' -insphere ' diam...
          ' -spec ' spec...
          ' -surf ' surf...
          ' -roi_nodes ' DpDa_math_60 '_node.txt'... 
          ' -roi_labels ' DpDa_math_60 '_label.txt'...
          ' -prefix std.60.rh.PP19_Dp-Da_math.MNI152.votc.inflated.' diam 'mm_diam']);
      
    spec = ['std.141.' sid '_rh.spec']; % Surface spec file
    surf = 'std.141.rh.inflated.gii'; % Surface to "grow" ROIs on
    unix(['ROIgrow -overwrite'... 
          ' -insphere ' diam...
          ' -spec ' spec...
          ' -surf ' surf...
          ' -roi_nodes ' DpDa_math_141 '_node.txt'... 
          ' -roi_labels ' DpDa_math_141 '_label.txt'...
          ' -prefix std.141.rh.PP19_Dp-Da_math.MNI152.votc.inflated.' diam 'mm_diam']);
end
  
