%% Convert Beta-series Pearson R-maps to p-value maps for conjunction analysis
purge
     
%% Get all R-maps
dir_start = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/';
cd(dir_start);
fnames = subdir('*Rmap*');

% Set the p-value thresholds for binary mask creation
p_thresholds = [0.05,0.01,0.005];

%% Loop through files
parfor ii = 1:numel(fnames)
    % Update on progress
    disp(['Working on file ' num2str(ii) ' of ' num2str(numel(fnames))]);
    
    % Read in the data
    f = fnames(ii).name;
    fdir = fnames(ii).folder;
    % Read in the events counts (number of valid trials for each condition)
    e = dir([fdir '/*_event_counts.txt']);
    ecounts = readmatrix([e.folder '/' e.name]);
    dset = afni_niml_readsimple(f);
    
    % Don't do the "OTHER" or "null" files for now
    if contains(f,'OTHER') || contains(f,'null')
        continue
    end
    
    % Set the number of events (N)
    % Check for the condition
    if contains(f,'Rmap.Da')
        N = ecounts(1,1);
    elseif contains(f,'Rmap.La')
        N = ecounts(2,1);
    elseif contains(f,'Rmap.Dp')
        N = ecounts(3,1);
    elseif contains(f,'Rmap.Lp')
        N = ecounts(4,1);
    elseif contains(f,'Rmap.ALL')
        N = sum(ecounts(1:4,1));
    end
    
    % Loop through all surface nodes and convert R to Pval
    dset_p = dset;
    for rr = 1:numel(dset.data)
        R = dset.data(rr);
        dset_p.data(rr) = convert_R_to_pval(R,N);
    end
    
    % Write out new dataset
    f_out = strrep(f,'Rmap','Pval');
    afni_niml_writesimple(dset_p,f_out);
    
    % Create thresholded binary maps
    for tt = 1:numel(p_thresholds)
        pt = p_thresholds(tt);
        pt_str = num2str(pt);
        dset_p_t = dset_p;
        dset_p_t.data = dset_p.data < pt;
        f_out = strrep(f,'Rmap',['Pval_thr' pt_str(3:end)]);
        afni_niml_writesimple(dset_p_t,f_out);
    end
end 

%% Global R-pval conversion
function pval = convert_R_to_pval(R,N)
    % R to P calculation
    t = sqrt(N-2)*R./sqrt(1-R.^2);
    s = tcdf(t,N-2);
    pval = 2 * min(s,1-s);
end











