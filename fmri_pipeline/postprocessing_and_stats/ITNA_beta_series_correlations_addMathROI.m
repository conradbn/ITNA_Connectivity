%% NFA Structure Function Project
% Beta-series correlations - Seed to wholebrain

purge
dir_start = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/';
cd(dir_start);
fnames = dir('*_proc');

%seed = {'Dp-Da' 'Lp-La'};
seed = {'Dp-Da_math'};
hemi = {'lh' 'rh'};

% For each subject
for sub = 1:numel(fnames)
    
    % Load the subject variables saved in workspace .mat file
    cd(dir_start)
    cd(fnames(sub).name);
    var_file = dir('workspace_variables_post*.mat');
    v = load(var_file.name,'subj','dir_results','dir_nii','dir_onset','dir_fs');
    
    % Correct folder paths to direct to server
    tmp = '/Users/nbl_imac/Documents/afni_tmp/';
    server = '/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/';
    dir_results = strrep(v.dir_results,tmp,server);
    dir_nii = strrep(v.dir_nii,tmp,server);
    dir_onset = strrep(v.dir_onset,tmp,server);
    dir_fs = strrep(v.dir_fs,tmp,server);
    subj = v.subj;
    
    % Set beta series results folder
    dir_beta = [dir_nii subj '.beta_series'];
    cd(dir_beta);
    
    %% Write event list to txt file
%     unix(['timing_tool.py -multi_timing ' dir_results 'stimuli/onsets_*_*_' subj '.txt'...
%         '  -multi_timing_to_event_list index ' subj '_event_list.txt'])
%     unix(['timing_tool.py -multi_timing ' dir_results 'stimuli/onsets_*_*_' subj '.txt'...
%         '  -multi_timing_to_event_list GE:ALL ' subj '_event_list_ALLinfo.txt'])
    
    % Load event list
    events = readmatrix([subj '_event_list.txt']);
    events = reshape(events',1,[])';
    events(isnan(events)) = [];
    
    % Write file with event counts
    ecounts = []; % because some subjects don't have error and/or ommission events, e.g. only 5 event types instead of 6
    [ecounts(:,1), ecounts(:,2)] = groupcounts(events);
%     writematrix(ecounts,[subj '_event_counts.txt']);
    
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
    
    % Set events to be censored (i.e., if more than 1 volume is censored
    % within the estimated HRF which is = 6 or 7 2s volumes)
    events_censor = num_TRs_censored > 1;
    
    %% Create beta-series correlation maps and z-score of difference
    label = {'Da' 'La' 'Dp' 'Lp'}; % Order is alphabetical based on onset filenames (Ca_Dig,Ca_Let,Cp_Dig,Cp_Let)
    contrast = {'Dp','Da';
                'Dp','Lp';
                'Lp','La'};
             
    % For seeds in each hemisphere
    for hh = 1:numel(hemi)
        h = hemi{hh};
        
        % Read in the beta series dataset
        betas = [dir_beta '/' subj '_' h '_LSS_betas.niml.dset'];
        b = afni_niml_readsimple(betas);
        
        % ******* MATH ROI ONLY IN RIGHT HEMISPHERE ********
        if strcmp(h,'lh')
            continue
        end
        % ******************************************
        
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
            for ee = 1:numel(label)
                e = label{ee};
                % For the data in each hemisphere
                for hhh = 1:numel(hemi)
                    hem = hemi{hhh};
                    %% True beta series correlations
                    % Read in the beta series dataset
                    betas = [dir_beta '/' subj '_' hem '_LSS_betas.niml.dset'];
                    betas_read = afni_niml_readsimple(betas);
                    % Get the index for this event type
                    idx = events == ee & events_censor == 0;
                    % Calculate correlation at every node
                    bs_corr = corr(bs_seed(idx),b.data(:,idx)');
                    % Write new niml surface dataset
                    bcorr = betas_read;
                    bcorr.data = bs_corr';
                    afni_niml_writesimple(bcorr,[subj '.' s '.' h '.beta_series_corr.' hem '.Rmap.' e '.niml.dset']);
                    % Fisher Z transform
                    bcorr.data = atanh(bcorr.data);
                    afni_niml_writesimple(bcorr,[subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' e '.niml.dset']);
                    
                    %% Null beta series correlations
                    % Now create a null dataset based on random set of
                    % betas for significance testing
                    % Get the index for this event type
                    events_all = find(events_censor ==0);
                    % Set up new dataset
                    bcorr_null = bcorr;
                    for nn = 1:100
                        idx_null = randperm(length(events_all),sum(idx));
                        idx_null = events_all(idx_null);
                        % Calculate correlation at every node
                        bs_corr = corr(bs_seed(idx_null),b.data(:,idx_null)');
                        bcorr_null.data(:,nn) = bs_corr';
                    end
                    afni_niml_writesimple(bcorr_null,[subj '.' s '.' h '.beta_series_corr.' hem '.Rmap.' e '_null.niml.dset']);
                    % Fisher Z transform
                    bcorr_null.data = atanh(bcorr_null.data);
                    afni_niml_writesimple(bcorr_null,[subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' e '_null.niml.dset']);
                     
                    %% Baseline beta series correlations 
                    % Now create mean connectivity map across ALL trials
                    % and all OTHER trials
                    % Only do the ALL trials map once
                    if ee == 1
                        % Get index for ALL trials (uncensored)
                        events_all = find(events_censor == 0);
                        % Calculate correlation at every node
                        bs_corr = corr(bs_seed(events_all),b.data(:,events_all)');
                        bcorr_all.data = bs_corr';
                        afni_niml_writesimple(bcorr_all,[subj '.' s '.' h '.beta_series_corr.' hem '.Rmap.ALL.niml.dset']);
                        % Fisher Z transform
                        bcorr_all.data = atanh(bcorr_all.data);
                        afni_niml_writesimple(bcorr_all,[subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.ALL.niml.dset']);
                    end
                    % Across OTHER trials
                    events_other = events ~= ee & events_censor == 0;
                    % Calculate correlation at every node
                    bs_corr = corr(bs_seed(events_other),b.data(:,events_other)');
                    bcorr_other.data = bs_corr';
                    afni_niml_writesimple(bcorr_other,[subj '.' s '.' h '.beta_series_corr.' hem '.Rmap.' e '_OTHER.niml.dset']);
                    % Fisher Z transform
                    bcorr_other.data = atanh(bcorr_other.data);
                    afni_niml_writesimple(bcorr_other,[subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' e '_OTHER.niml.dset']);
                    
                    %% Stimulus vs OTHER
                    % Get number of observations for each stimulus type
                    n1 = sum(events == ee & events_censor == 0);
                    n2 = sum(events_other);
                    % True data
                    z1 = afni_niml_readsimple([subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' e '.niml.dset']);
                    z2 = afni_niml_readsimple([subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' e '_OTHER.niml.dset']);
                    % Difference in Fisher Z values, taking into account degrees of freedom
                    zdiff = (z1.data - z2.data)./sqrt(1/(n1-3) + 1/(n2-3));
                    zout = z1;
                    zout.data = zdiff;
                    afni_niml_writesimple(zout,[subj '.' s '.' h '.beta_series_corr.' hem '.Zdiff.' e '-OTHER.niml.dset']);
                    
                end
            end
            
            % Create statistical contrast maps (Fisher-Z difference, as a z-score)
            % For each contrast
            for cc = 1:size(contrast,1)
                % Get condition labels for contrasting
                c1 = contrast{cc,1};
                c2 = contrast{cc,2};
                % Get number of observations for each stimulus type
                n1 = ecounts(strcmp(label,c1),1);
                n2 = ecounts(strcmp(label,c2),1);
                % For the data in each hemisphere
                for hhh = 1:numel(hemi)
                    hem = hemi{hhh};
                    % True data
                    z1 = afni_niml_readsimple([subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' c1 '.niml.dset']);
                    z2 = afni_niml_readsimple([subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' c2 '.niml.dset']);
                    % Difference in Fisher Z values, taking into account degrees of freedom
                    zdiff = (z1.data - z2.data)./sqrt(1/(n1-3) + 1/(n2-3));
                    zout = z1;
                    zout.data = zdiff;
                    afni_niml_writesimple(zout,[subj '.' s '.' h '.beta_series_corr.' hem '.Zdiff.' c1 '-' c2 '.niml.dset']);
                    
                    % Null data
                    z1 = afni_niml_readsimple([subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' c1 '_null.niml.dset']);
                    z2 = afni_niml_readsimple([subj '.' s '.' h '.beta_series_corr.' hem '.Zmap.' c2 '_null.niml.dset']);
                    % Difference in Fisher Z values, taking into account degrees of freedom
                    zdiff = (z1.data - z2.data)./sqrt(1/(n1-3) + 1/(n2-3));
                    zout = z1;
                    zout.data = zdiff;
                    afni_niml_writesimple(zout,[subj '.' s '.' h '.beta_series_corr.' hem '.Zdiff.' c1 '-' c2 '_null.niml.dset']);                
                end
            end
        end
    end
end