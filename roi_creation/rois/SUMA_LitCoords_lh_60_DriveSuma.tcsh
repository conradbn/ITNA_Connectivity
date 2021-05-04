#!/bin/tcsh -f
suma -spec /Users/benconrad/.afni/data/suma_MNI152_2009/std.60.MNI152_2009_lh.spec     -sv   /Users/benconrad/.afni/data/suma_MNI152_2009/MNI152_2009_SurfVol.nii -niml & sleep 5 & DriveSuma -com surf_cont -view_surf_cont y
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Abboud16_-40_-47_-12_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Carreiras15_-36_-62_-14_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Flowers04_-65_-60_-5_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Grotheer16_-47_-56_-14_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Pernet05_cat_-42.6_-61.8_-15.7_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Pernet05_disc_-41.5_-64.8_-13.8_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Polk02_actv_-45.75_-53_-6.75_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Polk02_pass_-41.7_-43.0_-8.7_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Pollack19_-42_-64_-11_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Puce96_mean_-40_-73.2_-19.2_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Rothlein14_-31_-59_-20_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Thesen12_-40_-78_-18_std.60_lh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
