%% Find closest surface nodes from MNI coordinates
% Using coordinates based on Pollack & Price 2019, NeuroImage
purge
MNI_dir = '/Users/benconrad/.afni/data/suma_MNI152_2009/';
cd(MNI_dir)
% Set the coordinates and labels
MNI_coord = {'-57 -52 -11'...
             '-42 -64 -11'...
             '57 -52 -11'...
             '42 -64 -11'};
label = {'Dp-Da_lh'...
         'Lp-La_lh'...
         'Dp-Da_rh'...
         'Lp-La_rh'};
% Loop through and write to text file
for ii = 1:numel(MNI_coord)
    m = MNI_coord{ii};
    l = label{ii};
    unix(['Surf2VolCoord -LPI'...
        ' -i ' MNI_dir 'std.60.lh.smoothwm.gii'...
        ' -i ' MNI_dir 'std.60.rh.smoothwm.gii'...
        ' -sv ' MNI_dir 'MNI152_2009_SurfVol.nii -qual LR'...
        ' -closest_node "' m '" > std.60.PP19_' l '.txt']);
end