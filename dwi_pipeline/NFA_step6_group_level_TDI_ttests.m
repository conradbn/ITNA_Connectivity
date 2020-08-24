purge
%% LH - Digit Region Projection Ends - from Pollack&Price19 ROIs
% Condition label
cmdlbl = 'Digit_TDI_ends_surf_lh_6mm';
% Outputs
topdir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
Outputs.dir = [topdir '/Group_StatisticalTests/'];
Outputs.maskgroup = [topdir '/Group_Masks/lh_rois.niml.dset[1]'];
Outputs.cmdlbl = cmdlbl;
Outputs.lbl =       {'Dp-Da'};
Outputs.script =    ['cmd_' cmdlbl '_3dttest++.txt'];
Outputs.output =    [cmdlbl '_3dttest++.nii.gz'];

% Inputs
Inputs.group_cmd = '3dttest++'; % '3dMEMA'
Inputs.betas =    '0';
Inputs.covars =   [topdir '/Covariates/covars_Nuisance.txt'];
Inputs.covarlbl = 'nuisance';
Inputs.sublist =  ''; %Empty uses all subjects
Inputs.blur = '0';%'0 2 6 6 8'; % Get outputs for each level of smoothing
Inputs.data1 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset']};     
Inputs.data2 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset']};

% RUN THE TTEST PIPELINE
NFA_group_level_TDI_ttest(Inputs,Outputs);

%% RH - Digit Region Projection Ends - from Pollack&Price19 ROIs
% Condition label
cmdlbl = 'Digit_TDI_ends_surf_rh_6mm';
% Outputs
topdir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
Outputs.dir = [topdir '/Group_StatisticalTests/'];
Outputs.maskgroup = [topdir '/Group_Masks/rh_rois.niml.dset[1]'];
Outputs.cmdlbl = cmdlbl;
Outputs.lbl =       {'Dp-Da'};
Outputs.script =    ['cmd_' cmdlbl '_3dttest++.txt'];
Outputs.output =    [cmdlbl '_3dttest++.nii.gz'];

% Inputs
Inputs.group_cmd = '3dttest++'; % '3dMEMA'
Inputs.betas =    '0';
Inputs.covars =   [topdir '/Covariates/covars_Nuisance.txt'];
Inputs.covarlbl = 'nuisance';
Inputs.sublist =  ''; %Empty uses all subjects
Inputs.blur = '0';%'0 2 6 6 8'; % Get outputs for each level of smoothing
Inputs.data1 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset']};     
Inputs.data2 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.rh.TDI_ends.norm.al2anat.lh.6mm.niml.dset']};

% RUN THE TTEST PIPELINE
NFA_group_level_TDI_ttest(Inputs,Outputs);

%% LH - Letter Region Projection Ends - from Pollack&Price19 ROIs
% Condition label
cmdlbl = 'Letter_TDI_ends_surf_lh_6mm';
% Outputs
topdir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
Outputs.dir = [topdir '/Group_StatisticalTests/'];
Outputs.maskgroup = [topdir '/Group_Masks/lh_rois.niml.dset[1]'];
Outputs.cmdlbl = cmdlbl;
Outputs.lbl =       {'Lp-La'};
Outputs.script =    ['cmd_' cmdlbl '_3dttest++.txt'];
Outputs.output =    [cmdlbl '_3dttest++.nii.gz'];

% Inputs
Inputs.group_cmd = '3dttest++'; % '3dMEMA'
Inputs.betas =    '0';
Inputs.covars =   [topdir '/Covariates/covars_Nuisance.txt'];
Inputs.covarlbl = 'nuisance';
Inputs.sublist =  ''; %Empty uses all subjects
Inputs.blur = '0';%'0 2 6 6 8'; % Get outputs for each level of smoothing
Inputs.data1 =     {[topdir '/1*/tracks_ss3t_*_combined_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.niml.dset']};     
Inputs.data2 =     {[topdir '/1*/tracks_ss3t_*_combined_Lp-La.lh.TDI_ends.norm.al2anat.rh.6mm.niml.dset']};

% RUN THE TTEST PIPELINE
NFA_group_level_TDI_ttest(Inputs,Outputs);

%% RH - Letter Region Projection Ends - from Pollack&Price19 ROIs
% Condition label
cmdlbl = 'Letter_TDI_ends_surf_rh_6mm';
% Outputs
topdir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
Outputs.dir = [topdir '/Group_StatisticalTests/'];
Outputs.maskgroup = [topdir '/Group_Masks/rh_rois.niml.dset[1]'];
Outputs.cmdlbl = cmdlbl;
Outputs.lbl =       {'Lp-La'};
Outputs.script =    ['cmd_' cmdlbl '_3dttest++.txt'];
Outputs.output =    [cmdlbl '_3dttest++.nii.gz'];

% Inputs
Inputs.group_cmd = '3dttest++'; % '3dMEMA'
Inputs.betas =    '0';
Inputs.covars =   [topdir '/Covariates/covars_Nuisance.txt'];
Inputs.covarlbl = 'nuisance';
Inputs.sublist =  ''; %Empty uses all subjects
Inputs.blur = '0';%'0 2 6 6 8'; % Get outputs for each level of smoothing
Inputs.data1 =     {[topdir '/1*/tracks_ss3t_*_combined_Lp-La.rh.TDI_ends.norm.al2anat.rh.6mm.niml.dset']};     
Inputs.data2 =     {[topdir '/1*/tracks_ss3t_*_combined_Lp-La.rh.TDI_ends.norm.al2anat.lh.6mm.niml.dset']};

% RUN THE TTEST PIPELINE
NFA_group_level_TDI_ttest(Inputs,Outputs);



    %% LH - Letter vs Letter Region Projection Ends Contrast - from Pollack&Price19 ROIs
    % Condition label
    cmdlbl = 'Digiti-Letter_TDI_ends_surf_lh_6mm_log';
    % Outputs
    topdir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
    Outputs.dir = [topdir '/Group_StatisticalTests/'];
    Outputs.maskgroup = [topdir '/Group_Masks/lh_rois.niml.dset[1]'];
    Outputs.cmdlbl = cmdlbl;
    Outputs.lbl =       {'Dp-Da';'Lp-La'};
    Outputs.script =    ['cmd_' cmdlbl '_3dttest++.txt'];
    Outputs.output =    [cmdlbl '_3dttest++.nii.gz'];

    % Inputs
    Inputs.group_cmd = '3dttest++'; % '3dMEMA'
    Inputs.betas =    '0';
    Inputs.covars =   [topdir '/Covariates/covars_Nuisance.txt'];
    Inputs.covarlbl = 'nuisance';
    Inputs.sublist =  ''; %Empty uses all subjects
    Inputs.blur = '0';%'0 2 6 6 8'; % Get outputs for each level of smoothing
    Inputs.data1 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log.niml.dset']
                       [topdir '/1*/tracks_ss3t_*_combined_Lp-La.lh.TDI_ends.norm.al2anat.lh.6mm.log.niml.dset']};     
    Inputs.data2 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.log.niml.dset']
                       [topdir '/1*/tracks_ss3t_*_combined_Lp-La.lh.TDI_ends.norm.al2anat.rh.6mm.log.niml.dset']};

    % RUN THE TTEST PIPELINE
    NFA_group_level_TDI_ttest(Inputs,Outputs);

    %% RH - Digit vs Letter Region Projection Ends Contrast - from Pollack&Price19 ROIs
    % Condition label
    cmdlbl = 'Digit-Letter_TDI_ends_surf_rh_6mm_log';
    % Outputs
    topdir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
    Outputs.dir = [topdir '/Group_StatisticalTests/'];
    Outputs.maskgroup = [topdir '/Group_Masks/rh_rois.niml.dset[1]'];
    Outputs.cmdlbl = cmdlbl;
    Outputs.lbl =       {'Dp-Da';'Lp-La'};
    Outputs.script =    ['cmd_' cmdlbl '_3dttest++.txt'];
    Outputs.output =    [cmdlbl '_3dttest++.nii.gz'];

    % Inputs
    Inputs.group_cmd = '3dttest++'; % '3dMEMA'
    Inputs.betas =    '0';
    Inputs.covars =   [topdir '/Covariates/covars_Nuisance.txt'];
    Inputs.covarlbl = 'nuisance';
    Inputs.sublist =  ''; %Empty uses all subjects
    Inputs.blur = '0';%'0 2 6 6 8'; % Get outputs for each level of smoothing
    Inputs.data1 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.rh.TDI_ends.norm.al2anat.rh.6mm.log.niml.dset']
                       [topdir '/1*/tracks_ss3t_*_combined_Lp-La.rh.TDI_ends.norm.al2anat.rh.6mm.log.niml.dset']};     
    Inputs.data2 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.rh.TDI_ends.norm.al2anat.lh.6mm.log.niml.dset']
                       [topdir '/1*/tracks_ss3t_*_combined_Lp-La.rh.TDI_ends.norm.al2anat.lh.6mm.log.niml.dset']};

    % RUN THE TTEST PIPELINE
    NFA_group_level_TDI_ttest(Inputs,Outputs);

        %% LH - Digit Projections vs Calc - from Pollack&Price19 ROIs
        covariates = {'Calc' 'Fluency' 'CalcSkills'};
        for ii = 1:numel(covariates)
            c = covariates{ii};
            % Condition label
            cmdlbl = ['Digit_x_' c '_TDI_ends_surf_lh_6mm_log'];
            % Outputs
            topdir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
            Outputs.dir = [topdir '/Group_StatisticalTests/'];
            Outputs.maskgroup = [topdir '/Group_Masks/lh_rois.niml.dset[1]'];
            Outputs.cmdlbl = cmdlbl;
            Outputs.lbl =       {'Dp-Da';'Lp-La'};
            Outputs.script =    ['cmd_' cmdlbl '_3dttest++.txt'];
            Outputs.output =    [cmdlbl '_3dttest++.nii.gz'];

            % Inputs
            Inputs.group_cmd = '3dttest++'; % '3dMEMA'
            Inputs.betas =    '0';
            Inputs.covars =   [topdir '/Covariates/covars_' c '.txt'];
            Inputs.covarlbl = 'nuisance';
            Inputs.sublist =  ''; %Empty uses all subjects
            Inputs.blur = '0';%'0 2 6 6 8'; % Get outputs for each level of smoothing
            Inputs.data1 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.lh.TDI_ends.norm.al2anat.lh.6mm.log.niml.dset']};     
            Inputs.data2 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.lh.TDI_ends.norm.al2anat.rh.6mm.log.niml.dset']};

            % RUN THE TTEST PIPELINE
            NFA_group_level_TDI_ttest(Inputs,Outputs);
        end

        %% RH - Digit Projections vs Calc - from Pollack&Price19 ROIs
        covariates = {'Calc' 'Fluency' 'CalcSkills'};
        for ii = 1:numel(covariates)
            c = covariates{ii};
            % Condition label
            cmdlbl = ['Digit_x_' c '_TDI_ends_surf_rh_6mm_log'];
            % Outputs
            topdir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152';
            Outputs.dir = [topdir '/Group_StatisticalTests/'];
            Outputs.maskgroup = [topdir '/Group_Masks/rh_rois.niml.dset[1]'];
            Outputs.cmdlbl = cmdlbl;
            Outputs.lbl =       {'Dp-Da';'Lp-La'};
            Outputs.script =    ['cmd_' cmdlbl '_3dttest++.txt'];
            Outputs.output =    [cmdlbl '_3dttest++.nii.gz'];

            % Inputs
            Inputs.group_cmd = '3dttest++'; % '3dMEMA'
            Inputs.betas =    '0';
            Inputs.covars =   [topdir '/Covariates/covars_' c '.txt'];
            Inputs.covarlbl = 'nuisance';
            Inputs.sublist =  ''; %Empty uses all subjects
            Inputs.blur = '0';%'0 2 6 6 8'; % Get outputs for each level of smoothing
            Inputs.data1 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.rh.TDI_ends.norm.al2anat.rh.6mm.log.niml.dset']};     
            Inputs.data2 =     {[topdir '/1*/tracks_ss3t_*_combined_Dp-Da.rh.TDI_ends.norm.al2anat.lh.6mm.log.niml.dset']};

            % RUN THE TTEST PIPELINE
            NFA_group_level_TDI_ttest(Inputs,Outputs);
        end
