function [] = NFA_group_level_TDI_ttest(Inputs,Outputs)

%% Set variables
% Start logging
diary([Outputs.dir 'log_' Outputs.cmdlbl '.txt']);
% Outputs
dir =       Outputs.dir;
cmdlbl =    Outputs.cmdlbl;
lbl =       Outputs.lbl;
script =    Outputs.script;
maskgroup = Outputs.maskgroup;

% Inputs
command =  Inputs.group_cmd;
betas =    Inputs.betas;
sublist =  Inputs.sublist;
covars =   Inputs.covars;
covarlbl = Inputs.covarlbl;
data1 =    Inputs.data1;
data2 =    Inputs.data2;

%% Set up output directory
cd(dir)
% Clear/Make output directory
unix(['rm -rf ' cmdlbl]);
mkdir(cmdlbl);  

% If covariates file, create string and copy to output folder
if isequal(covars,'')
    covars_str = '';
else 
    covars_str = [' -covariates ' covars];
    copyfile(covars,[dir '/' cmdlbl]);
end

% If subject list, add to command
if isequal(sublist,'')
    sublist_str = '';
else 
    sublist_str = [' -dset_sid_list ' sublist];
end

% Set output directory as the tempdir
tempdir = [dir cmdlbl];

% Copy mask to tempdir
unix(['cp ' maskgroup ' ' tempdir]);

%% Set hemisphere variables
if contains(cmdlbl,'lh')
    hemi = {'lh','rh'};
elseif contains(cmdlbl,'rh')
    hemi = {'rh','lh'};
end
%% Run tests for each hemisphere
for ii = 1:numel(hemi)
    h = hemi{ii};
    if ii == 1
        d = data1;
        mask_expr = 'a*not(b)'; % Use mask
    else
        d = data2;
        mask_expr = 'a'; % Dont use mask
    end
    
    % Generate group command for surface data
    cd(tempdir);
    s = [dir cmdlbl '/' h '_' script];
    if size(d,1) == 2
        unix(['gen_group_command.py -command ' command...
            ' -write_script ' s...
            ' -prefix ' cmdlbl '_3dttest++.' h '.niml.dset'...
            ' -dsets ' d{1}...
            ' ' sublist_str...
            ' -dsets ' d{2}...
            ' -set_labels ' lbl{1} ' ' lbl{2}...
            ' -subs_betas ' betas...
            ' -options \'...
            ' -paired ' covars_str...
            ' -tempdir ' tempdir]);
        %     ' -mask ' maskgroup ...
        %     ' -Clustsim '.
        %     ' -ETAC -ETAC_blur ' blur...
        %     ' -ETAC_opt NN=2:sid=2:hpow=0:pthr=0.005,0.001:name=ETAC_minimal']);

    elseif size(d,1) == 1
        unix(['gen_group_command.py -command ' command...
            ' -write_script ' s...
            ' -prefix ' cmdlbl '_3dttest++.' h '.niml.dset'...
            ' -dsets ' d{1}...
            ' ' sublist_str...
            ' -set_labels ' lbl{1} ...
            ' -subs_betas ' betas...
            ' -options \'...
            ' ' covars_str...
            ' -tempdir ' tempdir]);
        %     ' -Clustsim '...
        %     ' -mask ' maskgroup
        %     ' -ETAC -ETAC_blur ' blur...
        %     ' -ETAC_opt NN=2:sid=2:hpow=0:pthr=0.005,0.001:name=ETAC_minimal']);
    end

    % Run group level analysis command
    unix(['tcsh -ef ' s]);
    disp('++ Group level command executed!!');
    
    % Create masked result
    unix(['3dcalc -a ' cmdlbl '_3dttest++.' h '.niml.dset'...
        ' -b ' maskgroup '[1]'...
        ' -prefix ' cmdlbl '_3dttest++.' h '.mask.niml.dset'...
        ' -expr "' mask_expr '"']);
end

    
 
%% Generate cluster table for output
% unix(['3dClusterize -ithr 0 -idat 0 '...
%     ' -inset Ttest_Digit-Letter_TDI_ends_lh_3dttest++.default.ETACmaskALL.global.2sid.5perc.nii.gz[9]'...
%     ' -NN 2 -1sided RIGHT_TAIL 0.5 -clust_nvox 2 -pref_map ETAC.001.clust.order.nii.gz'])

%% Turn off logging and move log to output folder
diary off
movefile([Outputs.dir 'log_' Outputs.cmdlbl '.txt'],tempdir);

