clear all
% This script extracts the smoothing parameters from the initial runs of
% SurfSmooth (where smoothing was performed to a certain FWHM). I found
% that several subjects had many times too many iterations and resulted in
% extremely smooth data compared to the group. Instead I'll find the
% median parameters (iters and sigma) and hard code this into SurfSmooth
% for all data. 

cd('/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152')

% files = subdir('tracks*6mm.niml.dset.1D.smrec');
files = subdir('tracks*4mm.niml.dset.1D.smrec');

for ii = 1:numel(files) 
     fid = fopen(files(ii).name);
     cdata = textscan(fid,'%f%f','delimiter','\t', 'HeaderLines',6);
     fclose(fid);
     try
        outdata(ii,1) = cdata{1,1}(end);
     catch
        outdata(ii,1) = NaN;
     end
     try
        outdata(ii,2) = cdata{1,1}(end-1);
     catch
         outdata(ii,2) = NaN;
     end
     try
        outdata(ii,3) = cdata{1,2}(1);
     catch
         outdata(ii,3) = NaN;
     end
end
% Median sigma
median(outdata(:,1))
% Median Niters
median(outdata(:,2))

%% Save to file
%save('surface_smoothing_parameters_to_6mm_fwhm')
save('surface_smoothing_parameters_to_4mm_fwhm')

