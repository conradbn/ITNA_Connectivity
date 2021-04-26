clear all
cd('/Users/benconrad/Documents/GitHub/ITNA_Connectivity/roi_creation/rois')

dsets_lh = dir('LitCoord_*_std.141_lh.inflated.4mm_diam.niml.dset');
dsets_rh = dir('LitCoord_*_std.60_rh.inflated.4mm_diam.niml.dset');

%% LH
% mesh = '60';
% hemi = 'lh';
% dsets = dsets_lh;
% label = ['SUMA_LitCoords_' hemi '_' mesh];
% load_dsets_on_surface(mesh,hemi,label,dsets)

%% RH
mesh = '141';
hemi = 'rh';
dsets = dsets_rh;
label = 'LitCoord_RH';
load_dsets_on_surface(mesh,hemi,label,dsets)


%% Global function
function [] = load_dsets_on_surface(mesh,hemi,label,dsets)
    % Remove script if it exists
    unix(['rm -f ' label '_DriveSuma.tcsh']);
    % Set newline print command for cleaner looking code
    addstr = "fprintf(fid, '%s\n', str);";

    % Start new script file
    fid = fopen([label '_DriveSuma.tcsh'],'w');
    % Make script executable
    str = '#!/bin/tcsh -f'; eval(addstr);
    str = ['suma -spec /Users/benconrad/.afni/data/suma_MNI152_2009/std.' mesh '.MNI152_2009_' hemi '.spec'...
        '     -sv   /Users/benconrad/.afni/data/suma_MNI152_2009/MNI152_2009_SurfVol.nii -niml &'...
        ' sleep 5 &'... % insert a pause to allow SUMA to fully open
        ' DriveSuma -com surf_cont -view_surf_cont y'];eval(addstr);

    for ii = 1:numel(dsets)
        n = dsets(ii).name;
        if contains(n,'Letter')
            % Load the statistic dataset
            str = ['DriveSuma -com surf_cont -load_dset ' n ...
                ' -switch_cmap blue_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n'];eval(addstr);
        elseif contains(n,'Word')
            % Load the statistic dataset
            str = ['DriveSuma -com surf_cont -load_dset ' n ...
                ' -switch_cmap green_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n'];eval(addstr);
        elseif contains(n,'Digit')
            % Load the statistic dataset
            str = ['DriveSuma -com surf_cont -load_dset ' n ...
                ' -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n'];eval(addstr);
        end
    end

    % Kill the SUMA window
    %str = 'DriveSuma -com kill_suma & sleep 3'; eval(addstr);
    fclose(fid);

    % Actually run the final script
    %unix(['tcsh ' label '_DriveSuma.tcsh']);
end