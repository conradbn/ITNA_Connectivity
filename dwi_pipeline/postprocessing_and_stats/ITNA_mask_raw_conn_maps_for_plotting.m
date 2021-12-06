

data_dir = 'NBL_Projects/Price_NFA/Analyses_for_Paper/Results/Raw_Conn_Maps_Masked_for_Plotting';

% Group mean connectivity maps (after log transform)
conn_maps = {'PairedTest_SC_DigLH_vs_LetLH_log+c_SET1_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_LetLH_log+c_SET2_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_LetLH_log+c_contra_SET1_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_LetLH_log+c_contra_SET2_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf_SET1_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf_SET2_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf_contra_SET1_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf_contra_SET2_MEAN.niml.dset'};
         
% Consistency masks corresponding to each of the mean maps above         
masks = {'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'};


     
     