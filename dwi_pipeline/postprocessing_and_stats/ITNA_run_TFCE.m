%% Run Threshold-Free Cluster Enhancement (TFCE)
function [tfce_ds] = ITNA_run_TFCE(surf_ds,vertices,faces,niters,out,surf_ds_null)
    % All data is prepared; surf_ds has 8 samples and 5124 nodes. We want to
    % see if there are clusters that show a significant difference from zero in
    % their response. Thus, .sa.targets is set to all ones (the same
    % condition), whereas .sa.chunks is set to (1:8)', indicating that all
    % samples are assumed to be independent.
    %
    % (While this is a within-subject analysis, exactly the same logic can be
    % applied to a group-level analysis)

    % define neighborhood for each feature
    % (cosmo_cluster_neighborhood can be used also for meeg or volumetric
    % fmri datasets)
    cluster_nbrhood = cosmo_cluster_neighborhood(surf_ds,...
                                            'vertices',vertices,'faces',faces);

    fprintf('Cluster neighborhood:\n');
    cosmo_disp(cluster_nbrhood);

    opt=struct();

    % number of null iterations. for publication-quality, use >=1000;
    % 10000 is even better
    opt.niter = niters;

    % in this case we run a one-sample test against a mean of 0, and it is
    % necessary to specify the mean under the null hypothesis
    % (when testing classification accuracies, h0_mean should be set to chance
    % level, assuming a balanced crossvalidation scheme was used)
    
    if sum(unique(surf_ds.sa.targets)) == 1
        opt.h0_mean=0;
    end

    % this example uses the data itself (with resampling) to obtain cluster
    % statistcs under the null hypothesis. This is (in this case) somewhat
    % conservative due to how the resampling is performed.
    % Alternatively, and for better estimates (at the cost of computational
    % cost), one can generate a set of (say, 50) datasets using permuted data
    % e.g. using cosmo_randomize_targets), put them in a cell and provide
    % them as the null argument.
    if isempty(surf_ds_null)
        opt.null=[];
    else
        opt.null = surf_ds_null;
    end
    
    % Don't run stats testing if input is a single feature stat map
    %  'feature_stat' - 'none': Do not compute statistics, but instead
    %         use the input data from ds directly. In
    %         this case:
    %         * no statistic is computed;
    %         * ds.samples must be a row vector.
    %         * h0_mean is required.
    %         * the 'null' option is required
    %         * 'niter' must not be provided.
    %         * when using the 'tfce' cluster_stat
    %           option (the default), 'dh' must be
    %           provided explicitly.
    %         This option is intended for use when
    %         null data and feature-wise statistics
    %         have already been computed.    
    % 
    % 'dh',dh  (optional) Threshold step (only if cluster_stat is
    %       'tfce'). The default value of dh=0.1 should be fine
    %       in most (if not all) cases.
    %       Exception: when using 'feature_stat','none',
    %       the 'dh' option is required when using 'tfce'.
    %       For typical use cases, a value so that 100*dh is
    %       in the same order of magnitude as the range
    %       (maximum minus minimum) of the input (in .samples)
    %       may be a reasonable compromise between speed and
    %       accuracy
    
    if size(surf_ds.samples,1) == 1 && ~isempty(surf_ds_null)
        opt.feature_stat = 'none';
        opt.dh = 0.1;%0.02; %max(surf_ds.samples)-min(surf_ds.samples) ~= 10 ;
        opt = rmfield(opt,'niter');
    end
    
    
    % Use parallel computing
    opt.nproc = 8;
    fprintf('Running multiple-comparison correction with these options:\n');
    cosmo_disp(opt);
    
    % Run TFCE-based cluster correction for multiple comparisons.
    % The output has z-scores for each node indicating the probablity to find
    % the same, or higher, TFCE value under the null hypothesis
    tfce_ds=cosmo_montecarlo_cluster_stat(surf_ds,cluster_nbrhood,opt);


    %% Show results

    fprintf('TFCE z-score dataset\n');
    cosmo_disp(tfce_ds);

    nfeatures=size(tfce_ds.samples,2);
    percentiles=(1:nfeatures)/nfeatures*100;
    plot(percentiles,sort(tfce_ds.samples))
    title('sorted TFCE z-scores');
    xlabel('feature percentile');
    ylabel('z-score');


    nvertices=size(vertices,1);
    disp_opt=struct();
    disp_opt.DataRange=[-2 2];

    DispIVSurf(vertices,faces,1:nvertices,tfce_ds.samples',0,disp_opt);

    % % store results
    %fn_tfce_ds=fullfile(output_path, 'tfce_test1.niml.dset');
    cosmo_map2surface(tfce_ds,out);
    % 
    % surf_fn=fullfile(output_path, 'digit_intermediate.asc');
    % surfing_write(surf_fn,vertices,faces);
    % 
    % 
    % % show citation information
    % cosmo_check_external('-cite');
