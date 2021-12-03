%% ITNA Consistency Thresholding procedure
% Performs consistency-based thresholding to deal with false positives
% inherent in tractography results. Provides masks used for later
% statistical testing (i.e. only consider those nodes with reliable
% projections to/from seed)  Method inspired by Roberts et al. 2017
% (http://dx.doi.org/10.1016/j.neuroimage.2016.09.053)

purge;
cd('/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/AllSubs_dsets/')
% Set the files on which to run thresholding (this thresholding is
% performed on the raw TDI values before log transformation)
% fnames = {'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
%           'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'
%           'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm_MAP2CON.niml.dset'
%           'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
%           'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm_MAP2CON.niml.dset'
%           'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
%           'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'
%           'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.lh.6mm_MAP2CON.niml.dset'
%           'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.lh.6mm.niml.dset'};

% Input data, first lh, then rh
fnames = {'AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset','AllSubs_tracks_ss3t_50M_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset';
          'AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset','AllSubs_tracks_ss3t_50M_Lp-La.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset';
          'AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.lh.6mm.niml.dset','AllSubs_tracks_ss3t_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset'};
      
      
% Specify the ROIs to be used for masking      
% rois =  {'LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset'
%          'LitCoord_Letter_Pollack19_-42_-64_-11_std.141_lh.inflated.14mm_diam.niml.dset'
%          'LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam_MAP2CON.niml.dset'
%          'LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam.niml.dset'
%          'LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam_MAP2CON.niml.dset'};

rois =  {'LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset'
         'LitCoord_Letter_Pollack19_-42_-64_-11_std.141_lh.inflated.14mm_diam.niml.dset'
         'LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam.niml.dset'};
     
% rois_use_list = [1,2,3;...
%                  1,2,3;...
%                  1,2,3;...
%                  4,5,NaN;...
%                  4,5,NaN;...
%                  NaN,NaN,NaN;... % Contralateral hemisphere (no seed masking) 
%                  NaN,NaN,NaN;...
%                  NaN,NaN,NaN;...
%                  NaN,NaN,NaN];
rois_use_list = [1,NaN,NaN;...
                 2,NaN,NaN;...
                 3,NaN,NaN];
                
% Set proportion to keep based on consistency
prop_edges_keep_all = [0.005,0.01,0.02,0.04,0.08,0.15,0.30,0.60];% i.e. 0.30 = 30% density, based on edge consistency

% Set proportion of subjects required to have non-zero value.
% This is important because 0 values will be converted to NaN, and thus
% won't contribute the computation of std/mean. Even if there is a
% consistent edge in the remaining subjects, we don't want this to be in a
% small minority while 0 in rest of sample. (Note that conversion to NaN is
% important in case of log-transformed weights, which are negative).
% At least 90% of subjects need a nonzero value in this node (in this case, 26 of 29 subjects)
prop_subs_nzero = 0.90;

%% Loop through each of the files
for ii = 1:size(fnames,1)
    
    % Loop through each threshold level
    for pp = 1:numel(prop_edges_keep_all)
        prop_edges_keep = prop_edges_keep_all(pp);
        
        % Load the catenated subject data
        Ws_lh = afni_niml_readsimple(fnames{ii,1});
        Ws_rh = afni_niml_readsimple(fnames{ii,2});
        Ws = [Ws_lh.data;Ws_rh.data];
        
        % Mask the data with the ROIs so these nodes aren't considered in
        % consistency calculation (note these nodes will also be excluded from
        % later statistical testing)
        masks = rois_use_list(ii,:);
        for mm = 1:sum(~isnan(masks))
            r = rois{masks(mm)};
            m = afni_niml_readsimple(['/Users/nbl_imac2/Documents/GitHub/ITNA_Connectivity/roi_creation/rois/' r]);
            m = m.data(:,2);
            
            % Add zeroed data for the contralateral hemisphere
            if contains(r,'std.141_lh')
                m = [m; zeros(length(m),1)];
            elseif contains(r,'std.141_rh')
                m = [zeros(length(m),1); m];
            end
            % Multiply data by the inverse (binary) mask
            Ws = Ws .* ~m;
        end
        
        % Run the consistency thresholding process
        [W_thr,W_bin,W_log] = consistency_thresholding(Ws,prop_edges_keep,prop_subs_nzero);
        
        % Write thresholded and binarized surface files for each hemispher
        % seperately
        d_lh = Ws_lh;
        d_rh = Ws_rh;
        n_nodes = numel(d_lh.node_indices);
        
        d_lh.data = W_thr(1:n_nodes);
        afni_niml_writesimple(d_lh,['Consistency_Thresholding/consistency.combined_hemis.' num2str(prop_edges_keep) '.group_mean.' fnames{ii,1}]);
        d_lh.data = W_bin(1:n_nodes);
        afni_niml_writesimple(d_lh,['Consistency_Thresholding/consistency.combined_hemis.' num2str(prop_edges_keep) '.group_mask.' fnames{ii,1}]);
        
        d_rh.data = W_thr(n_nodes+1:end);
        afni_niml_writesimple(d_rh,['Consistency_Thresholding/consistency.combined_hemis.' num2str(prop_edges_keep) '.group_mean.' fnames{ii,2}]);
        d_rh.data = W_bin(n_nodes+1:end);
        afni_niml_writesimple(d_rh,['Consistency_Thresholding/consistency.combined_hemis.' num2str(prop_edges_keep) '.group_mask.' fnames{ii,2}]);
             
    end
    
    % Make sum images (i.e. sum of binary masks at each threshold)
    % LH
    all_lvls = dir(['Consistency_Thresholding/consistency.combined_hemis.0*group_mask*' fnames{ii,1}]);
    dall = [];
    for ll = 1:numel(all_lvls)
        d = afni_niml_readsimple(['Consistency_Thresholding/' all_lvls(ll).name]);
        dall(:,ll)=d.data;
    end
    d.data = nansum(dall,2);
    afni_niml_writesimple(d,['Consistency_Thresholding/consistency.combined_hemis.sum.group_mask.' fnames{ii,1}]);
    
    % RH 
    all_lvls = dir(['Consistency_Thresholding/consistency.combined_hemis.0*group_mask*' fnames{ii,2}]);
    dall = [];
    for ll = 1:numel(all_lvls)
        d = afni_niml_readsimple(['Consistency_Thresholding/' all_lvls(ll).name]);
        dall(:,ll)=d.data;
    end
    d.data = nansum(dall,2);
    afni_niml_writesimple(d,['Consistency_Thresholding/consistency.combined_hemis.sum.group_mask.' fnames{ii,2}]);
    
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
    W_bin = uint8(W_bin); % Set as integer data
    
    % Output a log10 version of the mean
    W_log = log10(Wmean);
    W_log(Wcv > t) = NaN;
end

