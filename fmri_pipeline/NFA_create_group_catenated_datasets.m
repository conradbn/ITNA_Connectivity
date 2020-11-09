%% Create group catenated/mean datasets to allow easy access for brain/behavior correlations
purge
out_dir = '/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations';

%% ROIs
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData';
cd(start_dir)

f = [start_dir '/*_proc/*.freesurfer/SUMA/std.141.lh.PollackandPrice19_Lp-La.MNI152.votc.inflated.14mm_diam.1.1D'];
make_mean_and_catenate_1D(f,out_dir)
f = [start_dir '/*_proc/*.freesurfer/SUMA/std.141.lh.PollackandPrice19_Dp-Da.MNI152.votc.inflated.14mm_diam.1.1D'];
make_mean_and_catenate_1D(f,out_dir)
f = [start_dir '/*_proc/*.freesurfer/SUMA/std.141.rh.PP19_Dp-Da_math.MNI152.votc.inflated.14mm_diam.1.1D'];
make_mean_and_catenate_1D(f,out_dir)
f = [start_dir '/*_proc/*.freesurfer/SUMA/std.60.lh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.1.1D'];
make_mean_and_catenate_1D(f,out_dir)
f = [start_dir '/*_proc/*.freesurfer/SUMA/std.60.lh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam.1.1D'];
make_mean_and_catenate_1D(f,out_dir)
f = [start_dir '/*_proc/*.freesurfer/SUMA/std.60.rh.PP19_Dp-Da_math.MNI152.votc.inflated.14mm_diam.1.1D'];
make_mean_and_catenate_1D(f,out_dir)

%% Activation
% Get the activation coefficients and tstats datasets
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData';
hemi = {'lh' 'rh'};
for ii = 1:numel(hemi)
    h = hemi{ii};
    f = [start_dir '/*_proc/*.results/stats.1*.' h '_REML.niml.dset'];
    make_mean_and_catenate(f,out_dir)
end

% Make datasets for each individual subbrik from full stats output
cd(out_dir)
task_conds = {'Dp' 'Da' 'Lp' 'La' 'Dp-Da' 'Lp-La'};

for ii = 1:numel(hemi)
    h = hemi{ii};
    % Load the stats datasets we just created
    stats = afni_niml_readsimple(['GroupMean_stats_' h '_REML.niml.dset']);
    all_stats = afni_niml_readsimple(['AllSubs_stats_' h '_REML.niml.dset']);
    for jj = 1:numel(task_conds)
        t = task_conds{jj};
        % Get the index for the Tstat and Beta (model coefficient) data
        brikT = find(strcmp(stats.labels,[t '#0_Tstat']));
        brikB = find(strcmp(stats.labels,[t '#0_Coef']));
        % Slice the data and write to files
        statsT = stats;
        statsT.data = all_stats.data(:,brikT:size(stats.data,2):end);
        statsB = stats;
        statsB.data = all_stats.data(:,brikB:size(stats.data,2):end);        
        afni_niml_writesimple(statsT,['AllSubs_stats_' h '_REML_Tstat_' t '.niml.dset']);
        afni_niml_writesimple(statsB,['AllSubs_stats_' h '_REML_Coeff_' t '.niml.dset']);
        % Calculate the the group mean and write to files
        statsT.data = mean(statsT.data,2);
        statsB.data = mean(statsB.data,2);
        afni_niml_writesimple(statsT,['GroupMean_stats_' h '_REML_Tstat_' t '.niml.dset']);
        afni_niml_writesimple(statsT,['GroupMean_stats_' h '_REML_Coeff_' t '.niml.dset']);
    end
end

%% Functional Connectivity
start_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData';
cd(start_dir)
hemi = {'lh' 'rh'};  
% Single conditions
task_conds = {'Dp' 'Da' 'Lp' 'La'};
for ii = 1:numel(hemi)
    for jj = 1:numel(task_conds)
        t = task_conds{jj};
        h = hemi{ii};
        f = [start_dir '/*_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' h '.Zmap.' t '.niml.dset'];
        make_mean_and_catenate(f,out_dir)
        f = [start_dir '/*_proc/*.beta_series/*.Lp-La.lh.beta_series_corr.' h '.Zmap.' t '.niml.dset'];
        make_mean_and_catenate(f,out_dir)
        f = [start_dir '/*_proc/*.beta_series/*.Dp-Da_math.rh.beta_series_corr.' h '.Zmap.' t '.niml.dset'];
        make_mean_and_catenate(f,out_dir)
    end
end

% Condition contrasts
task_conds = {'Dp-Da' 'Lp-La'};
for ii = 1:numel(hemi)
    for jj = 1:numel(task_conds)
        t = task_conds{jj};
        h = hemi{ii};
        f = [start_dir '/*_proc/*.beta_series/*.Dp-Da.lh.beta_series_corr.' h '.Zdiff.' t '.niml.dset'];
        make_mean_and_catenate(f,out_dir)
        f = [start_dir '/*_proc/*.beta_series/*.Lp-La.lh.beta_series_corr.' h '.Zdiff.' t '.niml.dset'];
        make_mean_and_catenate(f,out_dir)
        f = [start_dir '/*_proc/*.beta_series/*.Dp-Da_math.rh.beta_series_corr.' h '.Zdiff.' t '.niml.dset'];
        make_mean_and_catenate(f,out_dir)
    end
end

%% Structural connectivity
start_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface';
cd(start_dir)
hemi = {'lh' 'rh'};  
for ii = 1:numel(hemi)
    h = hemi{ii};
    f = ['*/tracks_ss3t_*_50M.wholebrain_TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M.wholebrain_TDI_ends.norm.al2anat.' h '.6mm.log.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M.wholebrain_length_map.al2anat.' h '.6mm.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M_Dp-Da.lh.TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M_Lp-La.lh.TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M_Dp-Da.lh.TDI_ends.norm.al2anat.' h '.6mm.log+c.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M_Lp-La.lh.TDI_ends.norm.al2anat.' h '.6mm.log+c.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M_Dp-Da_math.rh.TDI_ends.norm.al2anat.' h '.6mm.log+c.niml.dset'];
    make_mean_and_catenate(f,out_dir)
end

%% Global Functions
function make_mean_and_catenate_1D(f,out_dir)
    % First convert all the 1D files to niml
    files = dir(f);
    for ff = 1:numel(files)
        in = files(ff).name;
        in_dir = files(ff).folder;
        output = strrep(files(ff).name,'.1.1D','');
        cd(in_dir)
        if contains(f,'std.60.')
            unix(['ConvertDset -o_niml -pad_to_node ld60'...
                ' -input ' in ' -prefix ' output ' -node_index_1D ' in '[0]']);
        else
            unix(['ConvertDset -o_niml -pad_to_node ld141'...
                ' -input ' in ' -prefix ' output ' -node_index_1D ' in '[0]']);
        end
    end
    f = strrep(f,'.1.1D','.niml.dset');
    out = strsplit(f,'/');
    out = out{end};
    out = strrep(out,'_*_','_');
    out = strrep(out,'.1*.','_');
    out = strrep(out,'*.','');
    unix(['3dMean -non_zero -overwrite -prefix ' out_dir '/GroupMean_' out ' ' f]);
    unix(['3dTcat -overwrite -prefix ' out_dir '/AllSubs_' out ' ' f]);
end

function make_mean_and_catenate(f,out_dir)
    out = strsplit(f,'/');
    out = out{end};
    out = strrep(out,'_*_','_');
    out = strrep(out,'.1*.','_');
    out = strrep(out,'*.','');
    unix(['3dMean -non_zero -overwrite -prefix ' out_dir '/GroupMean_' out ' ' f]);
    unix(['3dTcat -overwrite -prefix ' out_dir '/AllSubs_' out ' ' f]);
end

