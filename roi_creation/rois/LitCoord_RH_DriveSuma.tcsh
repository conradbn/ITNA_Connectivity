#!/bin/tcsh -f
suma -spec /Users/benconrad/.afni/data/suma_MNI152_2009/std.60.MNI152_2009_rh.spec     -sv   /Users/benconrad/.afni/data/suma_MNI152_2009/MNI152_2009_SurfVol.nii -niml & sleep 5 & DriveSuma -com surf_cont -view_surf_cont y
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Gauthier00_54.8_-60.6_-3.7_std.60_rh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Longcamp04_40_-49_-14_std.60_rh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
DriveSuma -com surf_cont -load_dset LitCoord_Letter_Pernet05_disc_50.4_-66.2_-12.5_std.60_rh.inflated.4mm_diam.niml.dset -switch_cmap blue_monochrome -I_sb 1 -I_range 1 -Dim .6 -1_only n
