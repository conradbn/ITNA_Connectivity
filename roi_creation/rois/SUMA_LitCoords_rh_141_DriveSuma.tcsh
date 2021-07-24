#!/bin/tcsh -f
suma -spec /Users/benconrad/.afni/data/suma_MNI152_2009/std.141.MNI152_2009_rh.spec     -sv   /Users/benconrad/.afni/data/suma_MNI152_2009/MNI152_2009_SurfVol.nii -niml & sleep 5 & DriveSuma -com surf_cont -view_surf_cont y
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Abboud16_58_-46_-14_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Amalric16_62_-39_-17_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Bugden18_55_-51_-11_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Grotheer18_57_-54_-17_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Digit_Shum17_51_-54_-12_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap red_monochrome -I_sb 1 -I_range 0.5 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Gauthier00_54.8_-60.6_-3.7_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Longcamp04_40_-49_-14_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Pernet05_disc_50.4_-66.2_-12.5_std.141_rh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim 1 -1_only n
