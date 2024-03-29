#!/bin/tcsh -f
suma -spec /Users/benconrad/.afni/data/suma_MNI152_2009/std.141.MNI152_2009_lh.spec     -sv   /Users/benconrad/.afni/data/suma_MNI152_2009/MNI152_2009_SurfVol.nii -niml & sleep 5 & DriveSuma -com surf_cont -view_surf_cont y
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Amalric16_-56_-51_-19_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Bugden18_-53_-60_-9_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Grotheer18_-54_-55_-13_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Yeo17_-55_-50_-12_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Abboud16_-40_-47_-12_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Carreiras15_-36_-62_-14_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Flowers04_-65_-60_-5_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Grotheer16_-47_-56_-14_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Pernet05_cat_-42.6_-61.8_-15.7_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Pernet05_disc_-41.5_-64.8_-13.8_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Polk02_actv_-45.75_-53_-6.75_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Polk02_pass_-41.7_-43.0_-8.7_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Pollack19_-42_-64_-11_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Puce96_mean_-40_-73.2_-19.2_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Rothlein14_-31_-59_-20_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Thesen12_-40_-78_-18_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Word_Cattinelli13_-45_-47_-12_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap green_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Word_Cohen00_-42_-57_-6_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap green_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Word_Cohen04_-45_-57_-12_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap green_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Word_Dehaene10_-44_-50_-14_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap green_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Word_Thesen12_-46_-52_-20_std.141_lh.inflated.4mm_diam.niml.dset -switch_cmap green_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
