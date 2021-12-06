%% Combine ROIs and consistency thresholding masks for use in group-level tests
% Data within these masks would be excluded

purge
roi_dir = '/Users/nbl_imac2/Documents/GitHub/ITNA_Connectivity/roi_creation/rois';
consistency_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/AllSubs_dsets/Consistency_Thresholding';
cd(roi_dir)


% 'GroupMask_DigLH_LetLH_consistency_thresh_0.3_lh_ld141'
% 'GroupMask_DigLH_LetLH_consistency_thresh_0.3_rh_ld141'
% 'GroupMask_DigLH_DigRHflip_consistency_thresh_0.3_lh_ld141'
% 'GroupMask_DigLH_DigRHflip_consistency_thresh_0.3_rh_ld141'

%% Group masks for Digit vs Letter - Structural Connectivity 
% Includes consistency based threshold mask at 30% density (across both
% hemispheres combined) for each ROI. Either their combination (PLUS) or
% their conjunction (CONJ) Also masks the two ROIs themselves to avoid self
% connectivity in contrasts/stats.

% Left
label = 'GroupMask_DigLH_LetLH_consistency_thresh_0.3_lh';
m1 = afni_niml_readsimple([consistency_dir '/consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset']);
m2 = afni_niml_readsimple([consistency_dir '/consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset']);
r1 = afni_niml_readsimple([roi_dir '/LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset']);
r2 = afni_niml_readsimple([roi_dir '/LitCoord_Letter_Pollack19_-42_-64_-11_std.141_lh.inflated.14mm_diam.niml.dset']);

% Combined (PLUS) consistency mask and inverse
mout = m1;
mout.data = m1.data == 1 | m2.data == 1;
mout.data(r1.data(:,2) == 1 | r2.data(:,2) == 1) = 0;
afni_niml_writesimple(mout,[label '_PLUS.niml.dset']);

% Conjunction consistency mask
mout = m1;
mout.data = m1.data == 1 & m2.data == 1;
mout.data(r1.data(:,2) == 1 | r2.data(:,2) == 1) = 0;
afni_niml_writesimple(mout,[label '_CONJ.niml.dset']);


% Right (no ROIs in right hemisphere
label = 'GroupMask_DigLH_LetLH_consistency_thresh_0.3_rh';
m1 = afni_niml_readsimple([consistency_dir '/consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset']);
m2 = afni_niml_readsimple([consistency_dir '/consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset']);

% Combined (PLUS) consistency mask and inverse
mout = m1;
mout.data = m1.data == 1 | m2.data == 1;
afni_niml_writesimple(mout,[label '_PLUS.niml.dset']);

% Conjunction consistency mask
mout = m1;
mout.data = m1.data == 1 & m2.data == 1;
afni_niml_writesimple(mout,[label '_CONJ.niml.dset']);


%% Group masks for Digit Left vs Digit Right - Structural Connectivity 
% Includes consistency based threshold mask at 30% density (across both
% hemispheres combined) for each ROI. Either their combination (PLUS) or
% their conjunction (CONJ) Also masks the two ROIs themselves to avoid self
% connectivity in contrasts/stats.

% NOTE - the Digit RH connectivity data and masks have been mapped to the
% left hemispheres (MAP2CON)


% Left
label = 'GroupMask_DigLH_DigRH_consistency_thresh_0.3_lh';
m1 = afni_niml_readsimple([consistency_dir '/consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset']);
m2 = afni_niml_readsimple([consistency_dir '/consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm_MAP2CON.niml.dset']);
r1 = afni_niml_readsimple([roi_dir '/LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset']);
r2 = afni_niml_readsimple([roi_dir '/LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam_MAP2CON.niml.dset']);

% Combined (PLUS) consistency mask and inverse
mout = m1;
mout.data = m1.data == 1 | m2.data == 1;
mout.data(r1.data(:,2) == 1 | r2.data(:,2) == 1) = 0;
afni_niml_writesimple(mout,[label '_PLUS.niml.dset']);

% Conjunction consistency mask
mout = m1;
mout.data = m1.data == 1 & m2.data == 1;
mout.data(r1.data(:,2) == 1 | r2.data(:,2) == 1) = 0;
afni_niml_writesimple(mout,[label '_CONJ.niml.dset']);


% Right (no ROIs in right hemisphere
label = 'GroupMask_DigLH_DigRH_consistency_thresh_0.3_rh';
m1 = afni_niml_readsimple([consistency_dir '/consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset']);
m2 = afni_niml_readsimple([consistency_dir '/consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.lh.6mm_MAP2CON.niml.dset']);

% Combined (PLUS) consistency mask and inverse
mout = m1;
mout.data = m1.data == 1 | m2.data == 1;
afni_niml_writesimple(mout,[label '_PLUS.niml.dset']);

% Conjunction consistency mask
mout = m1;
mout.data = m1.data == 1 & m2.data == 1;
afni_niml_writesimple(mout,[label '_CONJ.niml.dset']);
      
   
   
   
   
