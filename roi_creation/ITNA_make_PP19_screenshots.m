%% Make PP19 ROI screenshots
% Use command line DriveSuma command to automatically start SUMA viewer
% load surface data, rotate as required, and generate image files for each
% subject

% * I have only been able to make this work in a Ubuntu 16.04 VM *

purge;

cd('/mnt/Price_NFA/NFA_fMRI/ProcessedData/');

subs = dir('1*');

for ii = 1:numel(subs)
    s = strsplit(subs(ii).name,'_');
    s = s{1};
    subFS_dir = ['/mnt/Price_NFA/NFA_fMRI/ProcessedData/' s '_proc/' s '.freesurfer/SUMA/'];
    cd(subFS_dir);
 
    
    d = 'std.141.lh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.1.1D';
    l = 'std.141.lh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam.1.1D';
    unix(['suma -spec std.141.13*_lh.spec -niml & '...
          ' DriveSuma '...
          ' -com surf_cont -load_dset ' d ' -switch_cmap red_monochrome -I_sb 1 -I_range 1'...
          ' -com surf_cont -load_dset ' l ' -switch_cmap blue_monochrome -I_sb 1 -I_range 1'...
          ' -com surf_cont -1_only n'... % Show both datasets
          ' -com viewer_cont -viewer_size 1000 1000'... % Window size
          ' -bkg_col 1 1 1'... % Background color (1 1 1 = white)
          ' -key:r2:s0.2 "."'... % Switch to inflated surface
          ' -key:s0.2 "a"'... % Make ROI's opaque
          ' -key:s0.2 "ctrl+left"'... % Go to left view
          ' -key:r7:s0.2 "up"'... % Slight rotation up to show lateral surface
          ' -key:r2:s0.2 "z"'... % Zoom out a bit
          ' -autorecord PP19_lh_MNI152_votc_coords_screenshots/PP19_lh_MNI152_votc_coords_dig-red_let-blue.inflated.14mm.' s... % Set image recording location
          ' -key:s0.2 "ctrl+r"']); % Actually save a screenshot, i.e. "record" 
    unix('DriveSuma -com kill_suma'); % Kill the SUMA window
    %' -key:r1:s0.2 "left"'... % Slight rotation left to show lateral surface
    cd('/mnt/Price_NFA/NFA_fMRI/ProcessedData/');
end


%% Copy all files to new folder
out_dir = '/mnt/Price_NFA/NFA_fMRI/ROI_Definition/PP19_screenshots';
fnames = subdir('PP19_lh_MNI152_votc_coords_dig-red_let-blue.inflated.14mm.*.jpg');
for ii = 1:numel(fnames)
    copyfile(fnames(ii).name,out_dir);
end

% Crop and flip the images, then catenate for each subject
cd(out_dir);
%unix('mogrify -flip *.jpg');
%unix('mogrify -flop *.jpg');
unix('mogrify -trim *MNI152*14mm*.jpg');


% ims = dir('lh*.jpg');
% for ii = 1:numel(ims)
%     s = strsplit(subs(ii).name,'.');
%     s = s{2};
%     %unix(['convert ' imgs(2).name ' ' imgs(1).name ' +append ' s '_combined.jpg']); 
%     %unix(['convert ' s '_combined.jpg -gravity southwest -pointsize 45 -annotate +0+0 "' s '" '  s '_combined_labeled.png']);
% end


unix('montage -mode concatenate -quality 100 PP19_lh_MNI152_votc_coords_dig-red_let-blue.inflated.14mm.*.jpg -tile 5x7 montage_PP19_lh_MNI152_votc_coords_dig-red_let-blue.inflated.14mm.png');
