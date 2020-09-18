
fnames = dir('1*/tracks_ss3t_*_50M_*p-*a.*h.TDI_ends.norm.al2anat.*h.6mm.niml.dset');
for ii = 1:numel(fnames)
    d = afni_niml_readsimple([fnames(ii).folder '/' fnames(ii).name]);
    d.data = log10(d.data + 0.5 * 1.4013e-45);
    afni_niml_writesimple(d,[fnames(ii).folder '/' fnames(ii).name(1:end-10) '.log+c.niml.dset']);
end


fnames = dir('1*/tracks_ss3t_*_50M_*p-*a.*h.TDI_ends.norm.al2anat.*h.6mm.log+c.niml.dset');
for ii = 1:numel(fnames)
    d = afni_niml_readsimple([fnames(ii).folder '/' fnames(ii).name]);
    dall(:,ii) = d.data;
end
