purge
fnames = dir('consistency.0*.group_mask.all_subs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');

for ii = 1:numel(fnames)
    d = afni_niml_readsimple(fnames(ii).name);
    dall(:,ii)=d.data;
end

d.data = nansum(dall,2);
afni_niml_writesimple(d,'consistency.sum.group_mask.all_subs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');




purge
fnames = dir('consistency.0*.group_mask.all_subs_tracks_ss3t_50M_Dp-Da.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset');

for ii = 1:numel(fnames)
    d = afni_niml_readsimple(fnames(ii).name);
    dall(:,ii)=d.data;
end

d.data = nansum(dall,2);
afni_niml_writesimple(d,'consistency.sum.group_mask.all_subs_tracks_ss3t_50M_Dp-Da.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset');



purge
fnames = dir('consistency.0*.group_mask.all_subs_tracks_ss3t_50M_Lp-La.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset');

for ii = 1:numel(fnames)
    d = afni_niml_readsimple(fnames(ii).name);
    dall(:,ii)=d.data;
end

d.data = nansum(dall,2);
afni_niml_writesimple(d,'consistency.sum.group_mask.all_subs_tracks_ss3t_50M_Lp-La.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset');




purge
fnames = dir('consistency.0*.group_mask.all_subs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');

for ii = 1:numel(fnames)
    d = afni_niml_readsimple(fnames(ii).name);
    dall(:,ii)=d.data;
end

d.data = nansum(dall,2);
afni_niml_writesimple(d,'consistency.sum.group_mask.all_subs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset');