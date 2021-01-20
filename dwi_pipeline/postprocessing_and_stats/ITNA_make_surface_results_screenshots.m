%% Make Screenshots of Final Connectivity Results
% Use command line DriveSuma command to automatically start SUMA viewer,
% load statistic datasets, load surface ROIs, rotate as required, and
% generate image files for each contrast

% NOTE - I was generating screenshots with DriveSuma in a Ubuntu 16.04 VM,
% but it was inconvenient.To make this finally work in OSX (High Sierra,
% but similar behavior on newer versions), I had to:
%
% 1) keep SUMA object controller permanently open as suggested here -
% https://afni.nimh.nih.gov/afni/community/board/read.php?1,155388,155419#msg-155419
%
% 2) resolve a dynamic library issue with the fix found here
% http://spedas.org/wiki/index.php?title=Known_Issues_macOS

%% Specify info about the contrast results of interest
purge;
top_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Results'; 
cd(top_dir);

input_strings = {
    % Structural Connectivity Contrasts
    'PairedTest_SC_DigLH_vs_LetLH_log+c','lh','ld141';...
    'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf','lh','ld141';...
    'PairedTest_SC_DigLH_vs_DigRH_log+c_on_RHsurf','rh','ld141';...

    % Functional Connectivity Contrasts - Digit vs Letter
    'PairedTest_FC_ALL_DigLH_vs_ALL_LetLH','lh','ld60';...
    'PairedTest_FC_Dp_DigLH_vs_Da_DigLH','lh','ld60';...
    'PairedTest_FC_Lp_LetLH_vs_La_LetLH','lh','ld60';...
    'PairedTest_FC_Dp-Da_DigLH_vs_Lp-La_LetLH','lh','ld60';...

    % Functional Connectivity Contrasts - Task-level        
    'PairedTest_FC_DTask_DigLH_vs_LTask_DigLH','lh','ld60';...
    'PairedTest_FC_DTask_LetLH_vs_LTask_LetLH','lh','ld60';...
    'PairedTest_FC_Dp_DigLH_vs_Lp_DigLH','lh','ld60';... 
    'PairedTest_FC_Dp_LetLH_vs_Lp_LetLH','lh','ld60';...
          
    % Functional Connectivity Contrasts - Digit L vs Digit R
    'PairedTest_FC_ALL_DigLH_vs_DigRH_on_LHsurf','lh','ld60';...
    'PairedTest_FC_ALL_DigLH_vs_DigRH_on_RHsurf','rh','ld60';...
    'PairedTest_FC_Dp_DigRH_vs_Da_DigRH_on_RHsurf','rh','ld60';...
    'PairedTest_FC_Dp-Da_DigLH_vs_DigRH_on_LHsurf','lh','ld60';...
    'PairedTest_FC_Dp-Da_DigLH_vs_DigRH_on_RHsurf','rh','ld60'};
   

%% Create screenshots and montage for each contrast result
for ii = 1:size(input_strings,1)
    % Get inputs 
    label = input_strings{ii,1};
    h = input_strings{ii,2};
    d = input_strings{ii,3};
    
    % Set strings
    dset = [label '_TFCE_Zscore_100000iters_MASK.niml.dset'];
    hemi = h;
    mesh = strrep(d,'ld','');
    cmap = 'coolwarm.niml.cmap';
    dimfac = '0.6';
    i_range = '4.26';
    t_thresh = '-T_val 1.96';
    
    % Set mask dataset(s) based on label
    dset_mask = {};
    if contains(label,'DigLH') && strcmp(hemi,'lh')
        dset_mask = [dset_mask,['LitCoord_Digit_Pollack19_-57_-52_-11_std.' mesh '_lh.inflated.14mm_diam.niml.dset']];
    elseif contains(label,'DigLH')
        dset_mask = [dset_mask,['LitCoord_Digit_Pollack19_-57_-52_-11_std.' mesh '_lh.inflated.14mm_diam_MAP2CON.niml.dset']];
    end
    if contains(label,'DigRH') && strcmp(hemi,'lh')
        dset_mask = [dset_mask,['LitCoord_Digit_Pollack19_54_-52_-14_std.' mesh '_rh.inflated.14mm_diam_MAP2CON.niml.dset']];
    elseif contains(label,'DigRH')
        dset_mask = [dset_mask,['LitCoord_Digit_Pollack19_54_-52_-14_std.' mesh '_rh.inflated.14mm_diam.niml.dset']];
    end
    if  contains(label,'LetLH') 
        dset_mask = [dset_mask,['LitCoord_Letter_Pollack19_-42_-64_-11_std.' mesh '_lh.inflated.14mm_diam.niml.dset']];
    end
    
    % Call functions to create screenshots (first remove those that exist)
    unix(['rm -f ' label '*.jpg & sleep 3']);
    call_SUMA(dset,label,hemi,mesh,cmap,dimfac,i_range,t_thresh,dset_mask);
    % Crop whitespace from all images and make consistent size
    unix(['mogrify -trim ' label '*.jpg & sleep 3']);
    unix(['mogrify -resize 1000x1000 ' label '*.jpg & sleep 3']);
    % Create montage of all four views
    create_montage(label,hemi);
end


%% Create screenshots and montages for raw connectivity maps
input_strings = {
    % Structural Connectivity
    'mask_PairedTest_SC_DigLH_vs_LetLH_log+c_SET1_MEAN.niml.dset','lh','ld141';...
    'mask_PairedTest_SC_DigLH_vs_LetLH_log+c_SET2_MEAN.niml.dset','lh','ld141';...
    'mask_PairedTest_SC_DigLH_vs_DigRH_log+c_on_RHsurf_SET2_MEAN.niml.dset','rh','ld141';...
    'PairedTest_SC_DigLH_vs_LetLH_log+c_SET1_MEAN.niml.dset','lh','ld141';...
    'PairedTest_SC_DigLH_vs_LetLH_log+c_SET2_MEAN.niml.dset','lh','ld141';...
    'PairedTest_SC_DigLH_vs_DigRH_log+c_on_RHsurf_SET2_MEAN.niml.dset','rh','ld141';...
    
    % Functional Connectivity
    'PairedTest_FC_ALL_DigLH_vs_ALL_LetLH_SET1_MEAN.niml.dset','lh','ld60';...
    'PairedTest_FC_ALL_DigLH_vs_ALL_LetLH_SET2_MEAN.niml.dset','lh','ld60';...
    'PairedTest_FC_ALL_DigLH_vs_DigRH_on_RHsurf_SET2_MEAN.niml.dset','rh','ld60'};

for ii = 4:6%1:size(input_strings,1)
    % Get inputs 
    dset = input_strings{ii,1};
    h = input_strings{ii,2};
    d = input_strings{ii,3};
    
    % Set strings
    label = strrep(dset,'.niml.dset','');
    hemi = h;
    mesh = strrep(d,'ld','');
    cmap = 'rainbow.niml.cmap';
    dimfac = '0.5';
    
    if contains(label,'SC')
        i_range = '-7 -4';
        t_thresh = '';
    elseif contains(label,'FC')
        i_range = '0 .5';
        t_thresh = '-T_val 0.141';
    end 
    
    % Set mask dataset(s) based on label
    if contains(label,'DigLH') && contains(label,'SET1')
       dset_mask = {['LitCoord_Digit_Pollack19_-57_-52_-11_std.' mesh '_lh.inflated.14mm_diam.niml.dset']};
    elseif contains(label,'LetLH') && contains(label,'SET2')
       dset_mask = {['LitCoord_Letter_Pollack19_-42_-64_-11_std.' mesh '_lh.inflated.14mm_diam.niml.dset']};
    elseif contains(label,'DigRH') && contains(label,'SET2')
       dset_mask = {['LitCoord_Digit_Pollack19_54_-52_-14_std.' mesh '_rh.inflated.14mm_diam.niml.dset']};
    end
    
    % Call functions to create screenshots (first remove those that exist)
    unix(['rm -f ' label '*.jpg & sleep 3']);
    call_SUMA(dset,label,hemi,mesh,cmap,dimfac,i_range,t_thresh,dset_mask);
    % Crop whitespace from all images and make consistent size
    unix(['mogrify -trim ' label '*.jpg & sleep 3']);
    unix(['mogrify -resize 1000x1000 ' label '*.jpg & sleep 3']);
    % Create montage of all four views
    create_montage(label,hemi);
end




%% Global function to call SUMA
% Will build a full tsch script line by line, then runs the sript
function call_SUMA(dset,label,hemi,mesh,cmap,dimfac, i_range,t_thresh,dset_mask) 
% Remove script if it exists
unix(['rm -f ' label '_DriveSuma.tcsh']);

% Set newline print command for cleaner looking code
addstr = "fprintf(fid, '%s\n', str)";

% Start new script file
fid = fopen([label '_DriveSuma.tcsh'],'w');

% Make script executable
str = '#!/bin/tcsh -f'; eval(addstr);

% Load MNI template surfaces  
% Also, keep SUMA object controller permanently open
% Necessary for DriveSuma to work in OSX, as suggested here:
% https://afni.nimh.nih.gov/afni/community/board/read.php?1,155388,155419#msg-155419
% MORE errors running this all from script - fix found here http://spedas.org/wiki/index.php?title=Known_Issues_macOS
str = ['suma -spec /Users/benconrad/.afni/data/suma_MNI152_2009/std.' mesh '.MNI152_2009_' hemi '.spec'...
      '     -sv   /Users/benconrad/.afni/data/suma_MNI152_2009/MNI152_2009_SurfVol.nii -niml &'...
      ' sleep 5 &'... % insert a pause to allow SUMA to fully open
      ' DriveSuma -com surf_cont -view_surf_cont y']; eval(addstr);

% Load the statistic dataset
str = ['DriveSuma -com surf_cont -load_dset ' dset ...
      ' -load_cmap ' cmap ' -I_sb 0 -I_range ' i_range...
      ' -T_sb 0 ' t_thresh ' -Dim ' dimfac ' -1_only n']; eval(addstr);

% Check if 1 or 2 masks specified and load
if numel(dset_mask) == 1
    str = ['DriveSuma -com surf_cont -load_dset ' dset_mask{1} ' -switch_cmap inverted_gray_circle -I_sb 1 -I_range 1']; eval(addstr);
elseif numel(dset_mask) == 2
    str = ['DriveSuma -com surf_cont -load_dset ' dset_mask{1} ' -switch_cmap inverted_gray_circle -I_sb 1 -I_range 1']; eval(addstr);
    str = ['DriveSuma -com surf_cont -load_dset ' dset_mask{2} ' -switch_cmap inverted_gray_circle -I_sb 1 -I_range 1']; eval(addstr);
end

% Inflate, rotate, etc, and record screenshots
str = ['DriveSuma -com viewer_cont -viewer_size 1000 1000'... % Window size
    ' -bkg_col 1 1 1'... % Background color
    ' -key:r2:s0.2 "."'... % Switch to inflated surface
    ' -key:s0.2 "a"'... % Make ROI's opaque
    ' -key:r2:s0.2 "z"'... % Zoom out a little bit
    ' -autorecord ' label ... % Set image recording name/path
    ' -key:s0.2 "ctrl+down" -key:s0.2 "ctrl+r"'... % Actually save a screenshots, i.e. "record"
    ' -key:s0.2 "ctrl+left" -key:s0.2 "ctrl+r"'...
    ' -key:s0.2 "ctrl+right" -key:s0.2 "ctrl+r"'...
    ' -key:s0.2 "ctrl+up" -key:s0.2 "ctrl+r"']; eval(addstr);

% Kill the SUMA window
str = 'DriveSuma -com kill_suma'; eval(addstr);
fclose(fid);

% Actually run the final script
unix(['./' label '_DriveSuma.tcsh']);
end

%% Function to create montage
function create_montage(label,hemi)
% Create montage of all views
imgs = dir([label '*.jpg']);
if strcmp(hemi,'lh')
    unix(['montage -tile 2x2 -mode Concatenate '...
        imgs(2).name ' ' imgs(3).name...
        ' \( ' imgs(4).name ' -rotate 90 \)'...
        ' \( ' imgs(1).name ' -rotate 270 \)'...
        ' ' label '_allviews_montage.jpg']);
elseif strcmp(hemi,'rh')
    unix(['montage -tile 2x2 -mode Concatenate '...
        imgs(3).name ' ' imgs(2).name...
        ' \( ' imgs(4).name ' -rotate 270 \)'...
        ' \( ' imgs(1).name ' -rotate 90 \)'...
        ' ' label '_allviews_montage.jpg']);
end
end


