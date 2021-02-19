%% Run Threshold Free Cluster Enhancement (TFCE) significance testing for all pairwise contrasts
purge
out_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Results';    

cd(out_dir)
niters = 100; % Total number of iterations for TFCE
test = 1; % All of the tests here are one-sample tests'

% input_strings = { % Linear Models
%                  'SC_DigL','lh','ld141',...
%                     'GroupMask_DigLH_LetLH_ConsThr0.3_ld141.niml.dset';            
%                  'SC_LetL','lh','ld141',...
%                     'GroupMask_DigLH_LetLH_ConsThr0.3_ld141.niml.dset' 
%                  'SC_DigR','rh','ld141',...
%                     'GroupMask_DigLH_CONJ_DigRH_ConsThr0.3_on_RHsurf_ld141.niml.dset'}; 
                 

input_strings = { % Linear Models
%                  'SC_DigL','lh','ld141',...
%                     'LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam_INV.niml.dset';
%                  'SC_DigR','rh','ld141',...
%                     'LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam_INV.niml.dset'
%                  'SC_LetL','lh','ld141',...
%                     'LitCoord_Letter_Pollack19_-42_-64_-11_std.141_lh.inflated.14mm_diam_INV.niml.dset'
                 'SC_DigL-DigR','lh','ld141',...
                    'LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam_INV.niml.dset';}; 
%                 'FC_DigL_ALL_Zmap','lh','ld60',...
%                     'LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset';      
%                  'FC_DigL_DpDa_Zdiff','lh','ld60',...
%                     'LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset';     
%                  'FC_DigL_LpLa_Zdiff','lh','ld60',...
%                     'LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam_INV.niml.dset';     
%                  'FC_DigR_ALL_Zmap','rh','ld60',...
%                     'LitCoord_Digit_Pollack19_54_-52_-14_std.60_rh.inflated.14mm_diam_INV.niml.dset';       
%                  'FC_DigR_DpDa_Zdiff','rh','ld60',...
%                     'LitCoord_Digit_Pollack19_54_-52_-14_std.60_rh.inflated.14mm_diam_INV.niml.dset';
                
dependent_vars = {'calc_skills_ss_resid_PP19'
                  'calc_skills_ss_PP19'
                  'calcss_resid_PP19'
                  'lwidss_resid_PP19'
                  'fluencyss_resid_PP19'
                  'ITG_CON_PP19'     % Note the ITG variables will be appended with L_ or R_ depending on the input data
                  'ITG_DBL_CON_PP19'};
              
                            
dependent_vars = {'LI_D_Pres'
                  'LI_CON_PP19'
                  'LI_D_Abs'};
              
%% 
for ii = 1:size(input_strings,1)
     cd(out_dir)
     data_dir = ['SemiPartCorr_' input_strings{ii,1}];
     hemi = input_strings{ii,2};
     density = input_strings{ii,3};
     mask_dset = input_strings{ii,4};
     cd(data_dir)
     
    for kk = 1:numel(dependent_vars)
        dv = dependent_vars{kk};
        % Set the hemisphere of the data
        if strcmp(hemi,'lh') && contains(dv,'ITG')
            dv = ['L_' dv];
        elseif strcmp(hemi,'rh') && contains(dv,'ITG')
            dv = ['R_' dv];
        end
 
        % Get the labels for this test       
        label = [out_dir '/' data_dir '/' strrep(data_dir,'SemiPartCorr_',['SemiPartCorr_' dv '_x_'])];
        
        % Single out the statistic dataset (i.e., the Tstat)
        dset = afni_niml_readsimple([label '.niml.dset']);
        dset.data = dset.data(:,2);
        in_dset = [label '_Tstat.niml.dset'];
        afni_niml_writesimple(dset,in_dset);
        
        % Set filename for the output statistic dataset and number of iterations
        out = [label '_TFCE_Zscore_' num2str(niters) 'iters_MASK.niml.dset'];
       
        % Specify the null data (i.e. the tstat datasets in which the dependent variable was shuffled)
        null_dsets = dir([label '*NULL_iter*.niml.dset']);
        null = dset;
        null.data = zeros(size(dset.data,1),numel(null_dsets));
        for jj = 1:numel(null_dsets)
            d = afni_niml_readsimple(null_dsets(jj).name);
            null.data(:,jj) = d.data(:,2); % 2nd column is the Tstat
        end
        out_dset_null = [label '_NULL_Tstats_ALL.niml.dset'];
        afni_niml_writesimple(null,out_dset_null);
        
        % Call global function
        setup_run_TFCE(in_dset,mask_dset,niters,out,test,hemi,density,out_dset_null)
      

           % Create null data based on sign-flipping (this is incorrect)
%         null = afni_niml_readsimple(in_dset);
%         truedata = null.data;
%         for jj = 1:niters
%             c = randi(2,size(truedata,1),1);
%             c(c == 2) = -1;
%             null.data(:,jj) = truedata.*c;
%         end
%         out_dset_null = [label '_NULL_Tstats_SIGNFLIP.niml.dset'];
%         afni_niml_writesimple(null,out_dset_null);
        
    end
end





%% Global function
function [] = setup_run_TFCE(in_dset,mask_dset,niters,out,test,hemi,density,out_dset_null)

% Load in the subject data
surf_ds = afni_niml_readsimple(in_dset);
% Mask the subject data
mask_ds = afni_niml_readsimple(['/Users/nbl_imac2/Documents/GitHub/ITNA_Connectivity/roi_creation/rois/' mask_dset]);
surf_ds.data = surf_ds.data .* mask_ds.data;

% Convert to COSMOMVPA formatted structure
surf_ds.labels = num2cell(1:size(surf_ds.data,2));
surf_ds = cosmo_surface_dataset(surf_ds);

% Get faces and vertices info from MNI surface gii file
if strcmp(density,'ld60')
    surf_gii = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/suma_MNI152_2009/std.60.' hemi '.smoothwm.gii'];
elseif strcmp(density,'ld141')
    surf_gii = ['/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/suma_MNI152_2009/std.141.' hemi '.smoothwm.gii'];
end

[vertices,faces]=surfing_read(surf_gii);
faces = double(faces);
vertices = double(vertices);

% Set up surface dataset structure, to prepare it for Cosmo TFCE
surf_ds.fa.center_ids = surf_ds.fa.node_indices;
if test == 2
    % Two-sample test
    nsubs = size(surf_ds.samples,1)/2;
    surf_ds.sa.chunks = [(1:nsubs)';(1:nsubs)'];
    surf_ds.sa.targets = [ones(nsubs,1);2*ones(nsubs,1)];
elseif test == 1
    % One-sample test (vs 0)
    nsubs = size(surf_ds.samples,1);
    surf_ds.sa.chunks = (1:nsubs)';
    surf_ds.sa.targets = ones(nsubs,1);
end

% Set up null surface data structure
surf_dset_null = cell(1,100);
if strcmp(out_dset_null,'') ~= 1
    null = afni_niml_readsimple(out_dset_null);
    ds = surf_ds;
    for ii = 1:size(null.data,2)%100
        %ds.samples = null.data(:,ii:100:size(null.data,2))';
        ds.samples = null.data(:,ii)';
        surf_dset_null(ii) = {ds};
    end
else
    surf_dset_null = [];
end

% Run TFCE
[~] = ITNA_run_TFCE(surf_ds,vertices,faces,niters,out,surf_dset_null);


% **** THE FOLLOWING ADDITIONAL TESTS DON'T MAKE SENSE IN THIS CONTEXT
% WHERE THE STATS DATASET IS THE INPUT AND WE ARE JUST APPLYING TFCE ******
%
% % Also get the uncorrected t-stats and z-scores for verification 
% % (NOTE this is adapted from last part of ITNA_run_TFCE)
% stat_ds = cosmo_stat(surf_ds,'t','z');
% tstat_ds = cosmo_stat(surf_ds,'t');
% fprintf('TFCE z-score dataset\n');
% cosmo_disp(stat_ds);
% nfeatures=size(stat_ds.samples,2);
% percentiles=(1:nfeatures)/nfeatures*100;
% plot(percentiles,sort(stat_ds.samples))
% title('sorted TFCE z-scores');
% xlabel('feature percentile');
% ylabel('z-score');
% nvertices=size(vertices,1);
% disp_opt=struct();
% disp_opt.DataRange=[-2 2];
% DispIVSurf(vertices,faces,1:nvertices,stat_ds.samples',0,disp_opt);
% 
% % Write uncorrected stats to surface files
% out_unc = strrep(out,'.niml.dset','_uncorrZscr.niml.dset');
% cosmo_map2surface(stat_ds,out_unc);
% out_uncT = strrep(out,'.niml.dset','_uncorrT.niml.dset');
% cosmo_map2surface(tstat_ds,out_uncT);

% Get Bayes Factors and write to surface file
% X = surf_ds.samples(1:nsubs,:);
% Y = surf_ds.samples(nsubs+1:end,:);
% bf10_ds = tstat_ds;
% parfor nn = 1:size(X,2)
%     x = X(:,nn)';
%     y = Y(:,nn)';
%     bf_xy(nn) = bf.ttest(x,y);   
% end
% bf10_ds.samples(1,:) = bf_xy;
% out_bf10 = strrep(out,'.niml.dset','_uncorrBF10.niml.dset');
% cosmo_map2surface(bf10_ds,out_bf10);

end
