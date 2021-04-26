#!/bin/tcsh -f
suma -spec /Users/benconrad/.afni/data/suma_MNI152_2009/std.60.MNI152_2009_rh.spec     -sv   /Users/benconrad/.afni/data/suma_MNI152_2009/MNI152_2009_SurfVol.nii -niml & sleep 5 & DriveSuma -com surf_cont -view_surf_cont y
