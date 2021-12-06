%% Mask Structural Connectivity Group Mean Maps
% Uses consistency thresholding maps at 30% density, plus absolute
% thresholds at log-7 and log-6, beyond which there is seemlingly little
% biological meaningfulness (cross-sectional area too small to consider).

clear all;
data_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Results/Raw_Conn_Maps_Masked_for_Plotting';
cd(data_dir)

% Group mean connectivity maps (after log transform)
conn_maps = {'PairedTest_SC_DigLH_vs_LetLH_log+c_SET1_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_LetLH_log+c_SET2_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_LetLH_log+c_contra_SET1_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_LetLH_log+c_contra_SET2_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf_SET1_MEAN.niml.dset'
             'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.log+c_MEAN.niml.dset'
             'PairedTest_SC_DigLH_vs_DigRH_log+c_on_LHsurf_contra_SET1_MEAN.niml.dset'
             'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.lh.6mm.log+c_MEAN.niml.dset'};
         
% Consistency masks corresponding to each of the mean maps above         
masks = {'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
         'consistency.combined_hemis.0.3.group_mask.AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'};
     
     
for ii = 1:numel(conn_maps)
    f = conn_maps{ii};
    m = afni_niml_readsimple(masks{ii});
    c = afni_niml_readsimple(f);
    
    masked = strrep(f,'.niml.dset','.Cons_MASK.niml.dset');
    masked_thresh7 = strrep(f,'.niml.dset','.Cons_MASK.grth-7.niml.dset');
    masked_thresh6 = strrep(f,'.niml.dset','.Cons_MASK.grth-6.niml.dset');
    
    % Consistency mask without absolute threshold
    cout = c;
    cout.data = cout.data .* m.data;
    afni_niml_writesimple(cout,masked);
    % Consistency mask + absolute threshold at -7
    cout.data(cout.data < -7) = NaN;
    afni_niml_writesimple(cout,masked_thresh7);
     % Consistency mask + absolute threshold at -6
    cout.data(cout.data < -6) = NaN;
    afni_niml_writesimple(cout,masked_thresh6);
    
end 



     
     