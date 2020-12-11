unix(['3dVol2Surf -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
' -sv ' suma_dir '/' sub '_SurfVol.nii'...
' -surf_A smoothwm'...
' -use_norms -norm_length -.5'...
' -f_index nodes'...
' -f_steps 10'...
' -map_func ave'...
' -oob_value 0'...
' -grid_parent ' tck '.TDI_ends.norm.al2anat.nii'...
' -out_niml h.-05.niml.dset']);
unix('3dhistog -omit 0 h.-05.niml.dset > h.-05.txt');

unix(['3dVol2Surf -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
' -sv ' suma_dir '/' sub '_SurfVol.nii'...
' -surf_A smoothwm'...
' -use_norms -norm_length -1'...
' -f_index nodes'...
' -f_steps 10'...
' -map_func ave'...
' -oob_value 0'...
' -grid_parent ' tck '.TDI_ends.norm.al2anat.nii'...
' -out_niml h.-10.niml.dset']);
unix('3dhistog -omit 0 h.-10.niml.dset > h.-10.txt');

unix(['3dVol2Surf -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
' -sv ' suma_dir '/' sub '_SurfVol.nii'...
' -surf_A smoothwm'...
' -use_norms -norm_length 1'...
' -f_index nodes'...
' -f_steps 10'...
' -map_func ave'...
' -oob_value 0'...
' -grid_parent ' tck '.TDI_ends.norm.al2anat.nii'...
' -out_niml h.10.niml.dset']);
unix('3dhistog -omit 0 h.10.niml.dset > h.10.txt');


unix(['3dVol2Surf -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
' -sv ' suma_dir '/' sub '_SurfVol.nii'...
' -surf_A smoothwm'...
' -use_norms -norm_length .5'...
' -f_index nodes'...
' -f_steps 10'...
' -map_func ave'...
' -oob_value 0'...
' -grid_parent ' tck '.TDI_ends.norm.al2anat.nii'...
' -out_niml h.05.niml.dset']);
unix('3dhistog -omit 0 h.05.niml.dset > h.05.txt');

unix(['3dVol2Surf -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
' -sv ' suma_dir '/' sub '_SurfVol.nii'...
' -surf_A smoothwm'...
' -f_index nodes'...
' -map_func mask'...
' -oob_value 0'...
' -grid_parent ' tck '.TDI_ends.norm.al2anat.nii'...
' -out_niml h.0.niml.dset']);
unix('3dhistog -omit 0 h.0.niml.dset > h.0.txt');

unix(['3dVol2Surf -spec ' suma_dir '/std.141.' sub '_' h '.spec'...
' -sv ' suma_dir '/' sub '_SurfVol.nii'...
' -surf_A smoothwm'...
' -surf_B pial'...
' -f_index nodes'...
' -f_steps 10'...
' -map_func ave'...
' -oob_value 0'...
' -grid_parent ' tck '.TDI_ends.norm.al2anat.nii'...
' -out_niml h.pial.niml.dset']);
unix('3dhistog -omit 0 h.pial.niml.dset > h.pial.txt');



unix('3dhistog -omit 0 tracks_ss3t_130843_combined_Dp-Da.lh.TDI_ends.norm.al2anat.mask2.nii.gz > h.vol.txt');


h1 = readtable('h.10.txt');
hn1 = readtable('h.-10.txt');
hn5 = readtable('h.-05.txt');
h5 = readtable('h.05.txt');
h0 = readtable('h.0.txt');
hpial = readtable('h.pial.txt');
hvol = readtable('h.vol.txt');

%% Plot
figure; hold on
plot(h1.x_Magnitude(2:35),h1.Freq(2:35))
plot(hn1.x_Magnitude(2:35),hn1.Freq(2:35))
plot(h5.x_Magnitude(2:35),h5.Freq(2:35))
plot(hn5.x_Magnitude(2:35),hn5.Freq(2:35))
plot(h0.x_Magnitude(2:35),h0.Freq(2:35))
plot(hpial.x_Magnitude(2:35),hpial.Freq(2:35))
plot(hvol.x_Magnitude(2:35),hvol.Freq(2:35))

legend('h1','hn1','h5','hn5','h0','hpial','hvol')

%% Average value
unix(['3dmaskave -quiet  h.10.niml.dset ']);
unix(['3dmaskave -quiet  h.05.niml.dset']);
unix(['3dmaskave -quiet   h.0.niml.dset ']);
unix(['3dmaskave -quiet   h.-05.niml.dset ']);
unix(['3dmaskave -quiet   h.-10.niml.dset  ']);
unix(['3dmaskave -quiet   h.pial.niml.dset ']);
unix(['3dmaskave -quiet   tracks_ss3t_130843_combined_Dp-Da.lh.TDI_ends.norm.al2anat.mask2.nii.gz']);


%% Average value
unix(['3dmaskave -quiet -dump  h.10.niml.dset > h10.txt']);
unix(['3dmaskave -quiet -dump  h.05.niml.dset > h05.txt']);
unix(['3dmaskave -quiet -dump  h.0.niml.dset > h0.txt']);
unix(['3dmaskave -quiet -dump  h.-05.niml.dset >  h-05.txt']);
unix(['3dmaskave -quiet -dump  h.-10.niml.dset > h-10.txt ']);
unix(['3dmaskave -quiet -dump  h.pial.niml.dset >  hpial.txt']);
unix(['3dmaskave -quiet -dump  tracks_ss3t_130843_combined_Dp-Da.lh.TDI_ends.norm.al2anat.mask2.nii.gz > hvol.txt']);

%% Plot histograms
close all
%t = load('hvol.txt'); t(t==0) = [];
x = load('h0.txt'); x(x==0) = [];
tt = hist(t,200);
xx = hist(x,200);
pdist2(tt(2:end),xx(2:end),'correlation');

