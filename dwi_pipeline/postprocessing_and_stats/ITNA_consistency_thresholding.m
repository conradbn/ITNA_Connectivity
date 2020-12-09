purge;
cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Data')
fnames = dir('all_subs*p-*a.*.TDI_ends.norm.al2anat.*.6mm.niml.dset');

% Set proportion to keep based on consistency
prop_edges_keep = 0.005 ; % i.e. 0.30 = 30% density, based on edge consistency

% Set proportion of subjects required to have non-zero value.
% This is important because 0 values will be converted to NaN, and thus
% won't contribute the computation of std/mean. Even if there is a
% consistent edge in the remaining subjects, we don't want this to be in a
% small minority while 0 in rest of sample. (Note that conversion to NaN is
% important in case of log-transformed weights, which are negative).
% At least 90% of subjects need a nonzero value in this node (in this case, 26 of 29 subjects)
prop_subs_nzero = 0.90;

%% Loop through each of the values
for ii = 1:numel(fnames)
    % Load the catenated subject data
    d = afni_niml_readsimple(fnames(ii).name);
    Ws = d.data;
    [W_thr,W_bin,W_log] = consistency_thresholding(Ws,prop_edges_keep,prop_subs_nzero);
    % Write thresholded and binarized surface files
    d.data = W_thr;
    afni_niml_writesimple(d,['consistency.' num2str(prop_edges_keep) '.group_mean.' fnames(ii).name]);
    d.data = W_bin;
    afni_niml_writesimple(d,['consistency.' num2str(prop_edges_keep) '.group_mask.' fnames(ii).name]);
    d.data = W_log.*W_bin;
    afni_niml_writesimple(d,['consistency.' num2str(prop_edges_keep) '.group_mean.log.' fnames(ii).name]);
end

%% Global function
function [W_thr,W_bin,W_log] = consistency_thresholding(Ws,prop_edges_keep,prop_subs_nzero)
    % Count how many subjects show a value of zero for each node
    node_zeros = sum(Ws == 0,2);
    nsubs = size(Ws,2);
    subthresh = nsubs*(1-prop_subs_nzero); 
    rm_ind = node_zeros >= subthresh;
    Ws(rm_ind,:) = 0;

    % Set all 0 values to NaN so they aren't included in computations
    Ws(Ws == 0) = NaN;

    % Compute mean weight matrix and CoV
    Wmean = nanmean(Ws,2); % group mean connectivity matrix
    Wcv = nanstd(Ws,0,2)./Wmean; % coefficient of variation across the group
    
%     if ~any(Wmean < 0)
%         Wcv = nanstd(Ws,0,2)./Wmean; % coefficient of variation across the group
%     end
%     else
%         Wcv = -nanstd(Ws,0,2)./Wmean; % coefficient of variation across the group
%     end

    % Find percentile (first convert NaN's to Inf because prctile will exclude
    % NaNs
    Wcv_inf = Wcv;
    Wcv_inf(isnan(Wcv_inf)) = Inf;
    %Wcv_inf(Wcv_inf == 0) = Inf;
    t = prctile(Wcv_inf,prop_edges_keep*100);

    % Create thresholded/binary weight vectors
    W_thr = Wmean;
    W_thr(Wcv > t) = NaN;
    W_bin = W_thr;
%     W_bin(W_bin ~= 0) = 1;
    W_bin(~isnan(W_bin)) = 1;
    
    % Output a log10 version of the mean
    W_log = log10(Wmean);
    W_log(Wcv > t) = NaN;
end

