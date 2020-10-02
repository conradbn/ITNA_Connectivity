%% NFA Structure Function Project
% Beta-series correlations - Seed to wholebrain

purge
dir_start = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/';
cd(dir_start);
fnames = dir('*_proc');

seed = {'Dp-Da' 'Lp-La'};
hemi = {'lh' 'rh'};

% For each subject
for sub = 1:numel(fnames)
    
    % Load the subject variables saved in workspace .mat file
    cd(dir_start)
    cd(fnames(sub).name);
    var_file = dir('workspace_variables_post*.mat');
    load(var_file.name,'subj','dir_results','dir_nii','dir_onset','dir_fs');
    
    % Correct folder paths to direct to server
    tmp = '/Users/nbl_imac/Documents/afni_tmp/';
    server = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/';
    dir_results = strrep(dir_results,tmp,server);
    dir_nii = strrep(dir_nii,tmp,server);
    dir_onset = strrep(dir_onset,tmp,server);
    dir_fs = strrep(dir_fs,tmp,server);
    
    % Set beta series results folder
    dir_beta = [dir_nii subj '.beta_series'];
    cd(dir_beta);
    
    %% Write event list to txt file
    unix(['timing_tool.py -multi_timing ' dir_results 'stimuli/onsets_*_*_' subj '.txt'...
        '  -multi_timing_to_event_list index ' subj '_event_list.txt'])
    unix(['timing_tool.py -multi_timing ' dir_results 'stimuli/onsets_*_*_' subj '.txt'...
        '  -multi_timing_to_event_list GE:ALL ' subj '_event_list_ALLinfo.txt'])
    
    % Load event list
    events = readmatrix([subj '_event_list.txt']);
    events = reshape(events',1,[])';
    events(isnan(events)) = [];
    
    %% Now, get the (per stimulus) censor data from afni files
    % Open afni_proc "review" txt output file to get number of censored TRs per stimulus trials
    fid  = fopen(['out.ss_review.' subj '_BSC.txt'],'r');
    text = textscan(fid,'%s','Delimiter','');
    text = text{1};
    fclose(fid);
    
    % Parse file for string of interest, i.e. number of TRs censored per
    % stimulus
    idx = contains(text,'num TRs censored per stim :');
    t = string(text(idx));
    t = regexp(t,'\d*','Match');
    num_TRs_censored = str2double(t');
   
    % Check that censor data and event list are same length
    if length(events) ~= length(num_TRs_censored)
        error('ERROR - Event list and censor info do not match in length');
    end
    
    % Set events to be censored
    events_censor = num_TRs_censored > 1;
    
    %% For each hemisphere
    for hh = 1:numel(hemi)
        h = hemi{hh};
        
        % Read in the beta series dataset
        betas = [dir_beta '/' subj '_' h '_LSS_betas.niml.dset'];
        b = afni_niml_readsimple(betas);
        
        % For each seed ROI
        for ss = 1:numel(seed)
            s = seed{ss};
            
            % Get the seed (ROI) mask node indices based on hemisphere/condition
            seed_mask = [dir_fs 'SUMA/std.60.' h '.PP19_' s '.MNI152.votc.inflated.14mm_diam.1.1D'];
            sm = readmatrix(seed_mask,'FileType','text','NumHeaderLines',2);
            
            % Get the beta series from the seed nodes
            bs = b.data(sm(:,1)+1,:); % Add one because AFNI index starts at 0
            bs_seed = mean(bs)';
            
            % Write to text file
            writematrix(bs',[s '.' h '.beta_series.all_nodes.1D'],'FileType','text');
            writematrix(bs_seed,[s '.' h '.beta_series.1D'],'FileType','text');
            
            % For the first 4 event types (Ca_Dig,Ca_Let,Cp_Dig,Cp_Let)
            label = {'Da' 'La' 'Dp' 'Lp'};
            for ee = 1:4
                % For the data in each hemisphere
                for hhh = 1:numel(hemi)
                    hem = hemi{hhh};
                    % Read in the beta series dataset
                    betas = [dir_beta '/' subj '_' hem '_LSS_betas.niml.dset'];
                    bb = afni_niml_readsimple(betas);
                    % Get the index for this event type
                    idx = events == ee & events_censor == 0;
                    % Calculate correlation at every node
                    bs_corr = corr(bs_seed(idx),b.data(:,idx)');
                    % Write new niml surface dataset
                    bcorr = bb;
                    bcorr.data = bs_corr';
                    afni_niml_writesimple(bcorr,[subj '.' s '.' h '.beta_series_corr.' hem '.Rmap.' label{ee} '.niml.dset']);
                    bcorr.data = atanh(bcorr.data);
                    afni_niml_writesimple(bcorr,[subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' label{ee} '.niml.dset']);
                end
            end
        end
    end
end