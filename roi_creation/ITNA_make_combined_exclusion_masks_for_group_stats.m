%% Combine ROIs and consistency thresholding masks for use in group-level tests
% Data within these masks would be excluded

purge
roi_dir = '/Users/nbl_imac2/Documents/GitHub/ITNA_Connectivity/roi_creation/rois';
cd(roi_dir)

%% First get combination and then inverse of consistency threshold masks
m1 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');
m2 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');
mout = m1;
mout.data = m1.data == 1 | m2.data == 1;
label = 'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_plus_LetLH.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset';
afni_niml_writesimple(mout,label);
mout.data = ~mout.data;
afni_niml_writesimple(mout,['INV_' label]);


m1 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');
m2 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm_MAP2CON.niml.dset');
mout = m1;
mout.data = m1.data == 1 | m2.data == 1;
label = 'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_plus_DigRH.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset';
afni_niml_writesimple(mout,label);
mout.data = ~mout.data;
afni_niml_writesimple(mout,['INV_' label]);


m1 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm_MAP2CON.niml.dset');
m2 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset');
mout = m1;
mout.data = m1.data == 1 | m2.data == 1;
label = 'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_plus_DigRH.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset';
afni_niml_writesimple(mout,label);
mout.data = ~mout.data;
afni_niml_writesimple(mout,['INV_' label]);

m1 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');
m2 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');
mout = m1;
mout.data = m1.data == 1 & m2.data == 1;
label = 'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_CONJ_LetLH.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset';
afni_niml_writesimple(mout,label);
mout.data = ~mout.data;
afni_niml_writesimple(mout,['INV_' label]);


m1 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');
m2 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm_MAP2CON.niml.dset');
mout = m1;
mout.data = m1.data == 1 & m2.data == 1;
label = 'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_CONJ_DigRH.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset';
afni_niml_writesimple(mout,label);
mout.data = ~mout.data;
afni_niml_writesimple(mout,['INV_' label]);


m1 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm_MAP2CON.niml.dset');
m2 = afni_niml_readsimple('consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset');
mout = m1;
mout.data = m1.data == 1 & m2.data == 1;
label = 'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_CONJ_DigRH.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset';
afni_niml_writesimple(mout,label);
mout.data = ~mout.data;
afni_niml_writesimple(mout,['INV_' label]);

%% Write the single ROI inverse maps (for masking the raw structural connectivity images)
label = {'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
        'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
        'consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'};
for ii = 1:numel(label)
    l = label{ii};
    m1 = afni_niml_readsimple(l);
    mout = m1;
    mout.data = mout.data~=1;
    afni_niml_writesimple(mout,['INV_' label{ii}]);
end


%% Now create combined masks for each density/contrast
% Functional connectivity masks
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_DigRH_on_LHsurf_ld60.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Digit_Pollack19_54_-52_-14_std.60_rh.inflated.14mm_diam_MAP2CON.niml.dset[1]'...
       ' -expr "not(or(a,b))"']);
   
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_DigRH_on_RHsurf_ld60.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_54_-52_-14_std.60_rh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam_MAP2CON.niml.dset[1]'...
       ' -expr "not(or(a,b))"']);  
   
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_LetLH_ld60.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Letter_Pollack19_-42_-64_-11_std.60_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -expr "not(or(a,b))"']);  

   
% Structural connectivity masks for contrasts   
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_DigRH_ConsThr0.3_on_LHsurf_ld141.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam_MAP2CON.niml.dset[1]'...
       ' -c INV_consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_plus_DigRH.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'...
       ' -expr "not(or(a,b,c))"']);
   
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_DigRH_ConsThr0.3_on_RHsurf_ld141.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam_MAP2CON.niml.dset[1]'...
       ' -c INV_consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_plus_DigRH.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'...
       ' -expr "not(or(a,b,c))"']);  
   
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_LetLH_ConsThr0.3_ld141.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Letter_Pollack19_-42_-64_-11_std.141_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -c INV_consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_plus_LetLH.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'...
       ' -expr "not(or(a,b,c))"']);  
   
   
% Structural connectivity masks for conjunctions    
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_CONJ_DigRH_ConsThr0.3_on_LHsurf_ld141.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam_MAP2CON.niml.dset[1]'...
       ' -c INV_consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_CONJ_DigRH.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'...
       ' -expr "not(or(a,b,c))"']);
   
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_CONJ_DigRH_ConsThr0.3_on_RHsurf_ld141.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam_MAP2CON.niml.dset[1]'...
       ' -c INV_consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_CONJ_DigRH.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'...
       ' -expr "not(or(a,b,c))"']);  
   
unix(['3dcalc -overwrite -prefix GroupMask_DigLH_CONJ_LetLH_ConsThr0.3_ld141.niml.dset'...
       ' -a LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -b LitCoord_Letter_Pollack19_-42_-64_-11_std.141_lh.inflated.14mm_diam.niml.dset[1]'...
       ' -c INV_consistency.0.3.group_mask.AllSubs_tracks_ss3t_50M_DigLH_CONJ_LetLH.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'...
       ' -expr "not(or(a,b,c))"']);  
      
      
   
   
   
   
