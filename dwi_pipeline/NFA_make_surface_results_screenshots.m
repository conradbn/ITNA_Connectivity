%% Make vOTC Searchlight results screenshots
% Use command line DriveSuma command to automatically start SUMA viewer
% load surface ROIs, rotate as required, and generate image files for each
% subject

% NOTE - This is not working properly in OSX. Seems that something is
% breaking with regards to interacting with SUMA from MATLAB. It is however
% working well in Ubuntu 16.04 VM.

purge;

top_dir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152/Group_StatisticalTests'; %'/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/'
cd(top_dir);
cons = dir('*_*');

%% Setup
dset = 'Digit-Letter_TDI_ends_surf_lh_6mm_log_3dttest++.lh.niml.dset';
label = 'Digit-Letter_TDI_ends_surf_lh_6mm_log';
hemi = 'lh';
mesh = '141';
cmap = 'Spectrum:yellow_to_cyan+gap';
i_range = '10';
t_thresh = '2.91';
clst_thresh = '601';
call_SUMA(dset,hemi,mesh,cmap,i_range,t_thresh)


% 
% %% Loop through each subject directory, call/drive SUMA, and save screenshots
% for ii = 1:numel(cons)
%     cd(top_dir);
%     c = cons(ii).name;
%     cd(c);
%     
%     if contains(c,'lh')
%         hemi = 'lh';
%     else
%         hemi = 'rh';
%     end
%     if contains(c,'smooth')
%         clst = '121';
%     else
%         clst = '10';
%     end
%     
%     % Right hemisphere
%     d = '*.niml.dset';
%     unix(['suma -spec MNI_N27:' hemi ':ld60 -niml & '...
%           ' DriveSuma '...
%           ' -com surf_cont -load_dset ' d ' -switch_cmap Spectrum:yellow_to_cyan+gap'...
%           ' -I_sb 1 -I_range 10 -T_sb 1 -T_val 2.91'... 
%           ' -com surf_cont -UseClst y -Clst -1 ' clst...'... % Show both datasets
%           ' -com viewer_cont -viewer_size 1000 1000'... % Window size
%           ' -bkg_col 1 1 1'... % Background color
%           ' -key:r2:s0.2 "."'... % Switch to inflated surface
%           ' -key:s0.2 "a"'... % Make ROI's opaque
%           ' -key:r2:s0.2 "z"'... % Zoom out a little bit
%           ' -autorecord ../screen_shot_' c ... % Set image recording location
%           ' -key:s0.2 "ctrl+down" -key:s0.2 "ctrl+r"'... % Actually save a screenshots, i.e. "record" 
%           ' -key:s0.2 "ctrl+left" -key:s0.2 "ctrl+r"'...
%           ' -key:s0.2 "ctrl+right" -key:s0.2 "ctrl+r"'...
%           ' -key:s0.2 "ctrl+up" -key:s0.2 "ctrl+r"']); 
%     unix('DriveSuma -com kill_suma'); % Kill the SUMA window
% %     ' -key:r3:s0.2 "z"'... % Zoom out a bit
% end
% 
% 
% %% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Copy all files to new folder
% out_dir = [top_dir '/GroupLevel/searchlight_screenshots'];
% fnames = subdir('*searchlight_dig-blue_let-red_postvOTC-green*.jpg');
% for ii = 1:numel(fnames)
%     copyfile(fnames(ii).name,out_dir);
% end
% 
% % Crop and flip the images, then catenate for each subject
% cd(out_dir);
% % unix('mogrify -flip *.jpg');
% % unix('mogrify -flop *.jpg');
% unix('mogrify -trim *.jpg');
% % for ii = 1:numel(cons)
% %     s = strsplit(cons(ii).name,'_');
% %     s = s{1};
% %     imgs = dir(['*' s '*']);
% %     unix(['convert ' imgs(2).name ' ' imgs(1).name ' +append ' s '_combined.jpg']); 
% %     unix(['convert ' s '_combined.jpg -gravity southwest -pointsize 45 -annotate +0+0 "' s '" '  s '_combined_labeled.png']);
% % end
% 
% unix('montage -mode concatenate -quality 100 *.jpg -tile 9x4 montage2.png');



%% Global function to call SUMA
function call_SUMA(dset,label,hemi,mesh,cmap,i_range,t_thresh)
    unix(['suma -spec /Users/benconrad/.afni/data/suma_MNI152_2009/std.' mesh '.MNI152_2009_' hemi '.spec'...
          ' -sv /Users/benconrad/.afni/data/suma_MNI152_2009/MNI152_2009_SurfVol.nii -niml &'...
          ' DriveSuma '...
          ' -com surf_cont -load_dset ' dset ' -switch_cmap ' cmap...
          ' -I_sb 1 -I_range ' i_range ' -T_sb 1 -T_val ' t_thresh... 
          ' -com surf_cont -UseClst y -Clst -1 ' clst_thresh...'... % Show both datasets
          ' -com viewer_cont -viewer_size 1000 1000'... % Window size
          ' -bkg_col 1 1 1'... % Background color
          ' -key:r2:s0.2 "."'... % Switch to inflated surface
          ' -key:s0.2 "a"'... % Make ROI's opaque
          ' -key:r2:s0.2 "z"'... % Zoom out a little bit
          ' -autorecord screen_shot_' label ... % Set image recording location
          ' -key:s0.2 "ctrl+down" -key:s0.2 "ctrl+r"'... % Actually save a screenshots, i.e. "record" 
          ' -key:s0.2 "ctrl+left" -key:s0.2 "ctrl+r"'...
          ' -key:s0.2 "ctrl+right" -key:s0.2 "ctrl+r"'...
          ' -key:s0.2 "ctrl+up" -key:s0.2 "ctrl+r"']); 
    unix('DriveSuma -com kill_suma'); % Kill the SUMA window
%     ' -key:r3:s0.2 "z"'... % Zoom out a bit
end



