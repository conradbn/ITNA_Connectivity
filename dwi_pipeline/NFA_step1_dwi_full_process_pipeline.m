function NFA_dwi_full_process_pipeline(sub_name,job_dir_path,group_proc_dir)
%% NFA DWI Full Process Pipeline (using SynB0-DISCO generated B0 image)
% This script runs the preprocessing for Price NFA DWI data using FSL's
% eddy tool via a pipeline by Justin Blaber of CCI/Masilab. Instead of B0
% fieldmap, here we use a synthetic B0 image which is corrected for
% distortion (using machine learning and the T1 image, via SynB0-Disco) as
% input into topup. The eddy outputs are then furthered processed to
% prepare for analysis in MRTrix3. 
% 10/2/2019

tic

%% Setup subject filepaths
%sub_name = 'Price_131031';
%job_dir_path = [sub_name '_eddy_no_field'];
raw_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/RawData/DICOM/';
proc_dir = [group_proc_dir job_dir_path];
field_dir = [proc_dir '/FIELDMAP'];

% Clear process folder 
unix(['rm -rf ' field_dir]);
unix(['rm -rf ' proc_dir]);

% Make subject process directory
unix(['mkdir ' proc_dir]);

% Write command window output to text file
diary([proc_dir '/command_window_output_' job_dir_path '.txt']);

%% Convert DICOMs to Nifti
cd(raw_dir);
fnames = dir([sub_name '*.DCM']);
for ii = 1:numel(fnames)
    dicm2nii(fnames(ii).name,field_dir,'.nii.gz');
end

% Go to fieldmap directory
cd(field_dir)

% Set input filenames
nii_anat = 'x3D_T1.nii.gz';
nii_dwi = 'WIP_HARDI_60_SENSE.nii.gz';
nii_fmap_real = 'Ax_B0map_real.nii.gz';
nii_fmap_mag = 'Ax_B0map_magnitude.nii.gz';

% Deoblique the raw nifti data (without this step, volumes are not aligned in
% viewers outside of AFNI, e.g. FSLeyes and MRIcron, not exactly sure
% why...). Affects only the header info, and does not touch the data. 
unix(['3drefit -deoblique ' nii_anat]);
unix(['3drefit -deoblique ' nii_dwi]);

%% Perform initial denoising of the DWI data (per MRTrix3 recommendation)
unix(['dwidenoise -force ' nii_dwi ' ' nii_dwi]);

%% Run SynB0 - Docker
% Set up input folder and files for SynB0 docker
synb0_dir = [proc_dir '/SYNB0'];
mkdir([synb0_dir '/INPUTS']);
unix(['3dcopy ' nii_anat ' ' synb0_dir '/INPUTS/T1.nii.gz']);
unix(['3dbucket -output ' synb0_dir '/INPUTS/b0.nii.gz ' nii_dwi '[0]']);
acqparams = [0,1,0,0.0144394;0,1,0,0.000];%0.062
writematrix(acqparams,[synb0_dir '/INPUTS/acqparams.txt'],'Delimiter',' ');

% Set FreeSurfer license path
fs_lic = '/Volumes/NBL_Projects/Price_NFA/Software+Code/dwi_pipeline/license.txt';

% Run SynB0 docker container
unix(['echo nbl_imac | sudo -S docker run --rm '...
     ' -v ' synb0_dir '/INPUTS/:/INPUTS/'...
     ' -v ' synb0_dir '/OUTPUTS/:/OUTPUTS/'...
     ' -v ' fs_lic ':/extra/freesurfer/license.txt'...
     ' --user $(id -u):$(id -g)'...
     ' justinblaber/synb0_25iso']);
 
%% Copy the SynB0 image to process folder 
copyfile([synb0_dir '/OUTPUTS/b0_u.nii.gz'],field_dir);
% Write bvec and bval file
bvec = [0;-0;0];
bval = 0;
dlmwrite([field_dir '/b0_u.bvec'],bvec);
dlmwrite([field_dir '/b0_u.bval'],bval);

%% EDDY PIPELINE ---------------------------------------------------------
% Run topup/eddy preprocessing pipeline written by Justin Blaber
% Use B0 fieldmap to correct distortions

% Set path to eddy pipeline utilities
addpath('/Volumes/NBL_Projects/Price_NFA/Software+Code/dwi_pipeline/system_utils');
addpath(genpath('/Volumes/NBL_Projects/Price_NFA/Software+Code/dwi_pipeline/nifti_utils'));
addpath(genpath('/Volumes/NBL_Projects/Price_NFA/Software+Code/dwi_pipeline/dwmri_visualizer'));
addpath('/Volumes/NBL_Projects/Price_NFA/Software+Code/dwi_pipeline/topup_eddy_preprocess');

% Start at group processing directory
cd(group_proc_dir);

% Set job directory path - ALREADY SET ABOVE
% job_dir_path = job_dir_path; %'nfa_test_topup_eddy_preprocess';

% Set FSL path
fsl_path = '/usr/local/fsl';

% BET params
bet_params = '-f 0.3 -R';

% Set dwmri_info - this will set base path to nifti/bvec/bval, phase 
% encoding direction, and readout times
dwmri_info(1).base_path = [field_dir '/' nii_dwi(1:end-7)];
dwmri_info(1).scan_descrip = 'scan';
dwmri_info(1).pe_dir = 'A';
dwmri_info(1).readout_time = 0.0144394;%0.062
dwmri_info(2).base_path = [field_dir '/b0_u'];
dwmri_info(2).scan_descrip = 'b0';
dwmri_info(2).pe_dir = 'P';
dwmri_info(2).readout_time = 0;

% ADC fix - apply it for Philips scanner
ADC_fix = true;

% zero_bval_thresh - will set small bvals to zero
zero_bval_thresh = 50;

% prenormalize - will prenormalize data prior to eddy
prenormalize = true;

% use all b0s for topup
use_all_b0s_topup = false;

% topup params
topup_params = ['--subsamp=1,1,1,1,1,1,1,1,1 ' ...
                '--miter=10,10,10,10,10,20,20,30,30 ' ...
                '--lambda=0.00033,0.000067,0.0000067,0.000001,0.00000033,0.000000033,0.0000000033,0.000000000033,0.00000000000067'];

% Sometimes name of eddy is 'eddy', 'eddy_openmp', or 'eddy_cuda'
eddy_name = 'eddy';%'eddy_openmp';

% use b0s in eddy
use_b0s_eddy = false;

% eddy params
eddy_params = '--repol --cnr_maps';
% if field_yes_no == 1
%     eddy_params = ['--repol --cnr_maps --field=' field_dir '/al_fieldmap_real_masked_reg'];
% elseif field_yes_no == 0
%     eddy_params = '--repol --cnr_maps';
% end

% normalize - will normalize data and output a single B0
normalize = true;

% sort scans - will sort scans by b-value
sort_scans = true;

% Set number of threads (only works if eddy is openmp version)
%setenv('OMP_NUM_THREADS','20');
setenv('OMP_NUM_THREADS','4');

% Perform preprocessing
[dwmri_path, bvec_path, bval_path, mask_path, movement_params_path, topup_eddy_pdf_path] = ...
    topup_eddy_preprocess(job_dir_path, ...
                          dwmri_info, ...
                          fsl_path, ...
                          ADC_fix, ...
                          zero_bval_thresh, ...
                          prenormalize, ...
                          use_all_b0s_topup, ...
                          topup_params, ...
                          eddy_name, ...
                          use_b0s_eddy, ...
                          eddy_params, ...
                          normalize, ...
                          sort_scans, ...
                          bet_params);
                      
%% Bias correction of DWI image
% *SKIP THIS STEP, DOING NORMALIZATION NOW WITH mtnormalize LATER IN PIPELINE*

% Bias correct the DWI images to remove B1 inhomogeneity This is
% recommended by MRTrix3 devs if ultimately running SIFT1/2 algorithms, as
% these would incorrectly interpret these inhomogenieties as
% spatially-varying fiber densities.

% % Go to eddy processed directory
% processed_dir = [proc_dir '/PREPROCESSED'];
% cd(processed_dir)
% 
% % Run bias correction
% nii_dwi_proc = 'dwmri.nii.gz'; 
% unix(['dwibiascorrect ' nii_dwi_proc ' bc_' nii_dwi_proc...
%      ' -ants -mask mask.nii.gz'...
%      ' -fslgrad dwmri.bvec dwmri.bval'...
%      ' -bias biasfield.nii.gz']);
% nii_dwi_proc = ['bc_' nii_dwi_proc];

%% Register T1/FreeSurfer to processed DWI
% First copy FreeSurfer data from SUMA folder (created by AFNI during fMRI
% surface analysis preprocessing)
sub_id = sub_name(7:end);
nii_anat = 'brain.finalsurfs.nii';
nii_seg = 'aseg.nii';
suma_dir = ['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sub_id '_proc/' sub_id '.freesurfer/SUMA/'];
copyfile([suma_dir nii_anat], processed_dir);
copyfile([suma_dir nii_seg], processed_dir);

% Get the dwi b0 image
unix(['3dbucket -prefix dwmri_b0.nii.gz ' nii_dwi_proc '[0]']);

%--------------------------------------------------------------------------
% LOAD THE ORIGINAL WORKSPACE
proc_dir = [group_proc_dir job_dir_path];
cd([proc_dir '/PREPROCESSED']);
% Files to keep
%keep = {'aseg.nii';'mask.nii.gz';'brain.finalsurfs.nii';'eddy_params.txt'}


% AFFINE ALIGNMENT --------------------------------------------------------
% Align the anatomical to dwi b0 image/space
unix(['align_epi_anat.py -suffix _al2dwi'...
      ' -epi_base 0 -anat2epi -ginormous_move -partial_axial'...
      ' -anat ' nii_anat ' -anat_has_skull no'...
      ' -master_epi ' nii_anat...
      ' -epi dwmri_b0.nii.gz -epi_strip 3dAutomask']);
% Apply transform to aseg data  
unix(['3dAllineate -final NN -1Dmatrix_apply brain.finalsurfs_al2dwi_mat.aff12.1D'...
      ' -master ' nii_seg...
      ' -input ' nii_seg...
      ' -prefix aseg_al2dwi.nii.gz']);
unix('3dcopy brain.finalsurfs_al2dwi+orig brain.finalsurfs_al2dwi.nii.gz');
% ------------------------------------------------------------------------- 
  
% AFFINE + NONLINEAR ALIGNMENT -------------------------------------------- 
% Align the dwi b0 to the anatomical image/space
unix(['align_epi_anat.py -suffix _al2anat'...
      ' -epi_base 0 -epi2anat -ginormous_move -partial_axial'... % -big_move (not enough to align some subjects)
      ' -anat ' nii_anat ' -anat_has_skull no'...
      ' -master_epi ' nii_anat...
      ' -epi dwmri_b0.nii.gz -epi_strip 3dAutomask']);

% Nonlinear alignment using only global warping. 
% DWI (EPI) image is used as the base, so it provides the weighting.
unix(['3dQwarp -source ' nii_anat...
      ' -base dwmri_b0_al2anat+orig'...
      ' -prefix anat_warp2tmp_dwi.nii.gz'...
      ' -lpc -verb -iwarp -blur 0 3']);

% Use the inverse warp to transform FreeSurfer data to DWI b0 space
unix(['3dNwarpApply -prefix aseg_nl2dwi.nii.gz'...
      ' -source ' nii_seg ' -ainterp NN'...
      ' -newgrid 1'...
      ' -nwarp "anat_warp2tmp_dwi_WARP.nii.gz brain.finalsurfs_al2anat_mat.aff12.1D"']);

unix(['3dNwarpApply -prefix brain.finalsurfs_nl2dwi.nii.gz'...
      ' -source ' nii_anat...
      ' -newgrid 1'...
      ' -nwarp "anat_warp2tmp_dwi_WARP.nii.gz brain.finalsurfs_al2anat_mat.aff12.1D"']);
% -------------------------------------------------------------------------
% Save workspace
save(job_dir_path);

%% MRTrix3 Processing 
% Generate the 5 tissue file for ACT (Anatomically constrained tractography) 
unix('5ttgen freesurfer -lut /Applications/freesurfer/FreeSurferColorLUT.txt aseg_al2dwi.nii.gz aseg_5ttgen.mif'); 

% Generate GM/WM interface for seeding
unix('5tt2gmwmi aseg_5ttgen.mif aseg_5tt_gmwmi.mif');

% Convert nifti/bvec/bval data to mif format (MRtrix image format)
unix(['mrconvert -fslgrad dwmri.bvec dwmri.bval ' nii_dwi_proc ' dwmri.mif']);

%% Turn off writing command window output to text file
toc
diary('off');





























 