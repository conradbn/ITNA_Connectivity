 
purge
cutoff = 1e-5;
a = afni_niml_readsimple('tracks_ss3t_132531_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset');
b = afni_niml_readsimple('tracks_ss3t_132531_50M_Lp-La.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset');

a(a.data) =
a = normalize(a.data);
b = normalize(b.data);

d = [a,b];
d(d<cutoff)=NaN;

figure; hold on

histogram(d(:,1))
histogram(d(:,2))




cosmo_check_external({'surfing','afni'});