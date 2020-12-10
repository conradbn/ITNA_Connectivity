%% Convert Beta-series Pearson R-maps to p-value maps for conjunction analysis
purge
     
%% Get all Zdiff-maps
dir_start = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/';
cd(dir_start);
fnames = subdir('*Zdiff*');

% Set the p-value thresholds for binary mask creation
p_thresholds = [0.05,0.01,0.005];


%% Loop through files
for ii = 1:numel(fnames)
    % Update on progress
    disp(['Working on file ' num2str(ii) ' of ' num2str(numel(fnames))]);
    
    % Read in the data
    f = fnames(ii).name;
    fdir = fnames(ii).folder;
    dset = afni_niml_readsimple(f);
    
    % Don't do the "OTHER" or "null" files for now
    if contains(f,'OTHER') || contains(f,'null')
        continue
    end
    
    % Loop through all surface nodes and convert Z-score of difference to Pval
    % DON'T do this because only care about positive effects for now
%     p = normcdf(dset.data);
%     p(p>=0.5) = 1-p;
%     dset_p.data = p;
    
    for tt = 1:numel(p_thresholds)
        dset_p = dset;
        thresh = norminv(1-p_thresholds(tt)); % One-sided test (i.e. for 0.05, use norminv(0.95))
        dset_p.data = dset.data > thresh;
        pt_str = num2str(p_thresholds(tt));
        f_out = strrep(f,'Zdiff',['Zdiff_Pval_thr' pt_str(3:end)]);
        afni_niml_writesimple(dset_p,f_out);
    end
end 










