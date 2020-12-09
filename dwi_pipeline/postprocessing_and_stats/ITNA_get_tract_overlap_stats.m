%% Get overlap of TDI and largescale tracts
% This script calculates overlap statistics of subject-level TDI maps (in
% MNI space), with all tracts/bundles defined in the (4D) Pandora White
% Matter Atlas @ https://www.nitrc.org/projects/pandora_atlas

purge 

atlas_names = {'AFQ' 'Xtract'};%'AFQclipped' 'Recobundles' 'TractSeg' 'Tracula' 
atlas_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/Pandora_Atlas';%'/Users/benconrad/Downloads/NITRC-multi-file-downloads'; 
mask = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_MNI152/MNI152_2009_template_mask.nii.gz';

% Get filenames for each subject
%cd /Users/nbl_imac2/Desktop/Price_NFA_Tractography_MNI152/
cd /Volumes/BensHD_2020/Price_NFA_Tractography_MNI152
%niis = subdir('Group_Average_*p-*a_rh_fwhm4.nii.gz'); %'tracks_ss3t_*_combined_Dp-Da_MNI_TDI_fwhm4.nii.gz');
niis = subdir('Digit-Letter_TDI_rh_4mm_3dttest++_Coef_masked_*.nii'); 
niis = subdir('tracks_ss3t_*_combined_Dp-Da_MNI_TDI_fwhm4.nii.gz');

%% Run across subjects/atlases
for ii = 1:numel(atlas_names)
    a_name = atlas_names{ii};
    a_table = readtable([atlas_dir '/' a_name '/' a_name '_info.csv']);
    a_nii = [atlas_dir '/' a_name '/' a_name '.nii.gz']; %load_untouch_nii_gz([atlas_dir '/' a_name '/' a_name '.nii.gz']);
    for jj = 1:numel(niis)
        n = niis(jj).name;
        for kk = 1:size(a_table,1)
            disp(['--------- WORKING ON '...
                  ' atlas ' num2str(ii) ' of ' num2str(numel(atlas_names)) ...
                  ' / subject ' num2str(jj) ' of ' num2str(numel(niis))...
                  ' / tract ' num2str(kk) ' of ' num2str(size(a_table,1))])
            % Setup inputs
            base = a_nii;
            base_ind = kk-1;
            source = n;
            source_mask = mask;
            %cost_fxn = 13; % LPC+, combination metric
            % Run function
            overlap.(a_name)(jj,kk,:) = calculate_overlap(base,base_ind,source,source_mask);
            overlap.([a_name '_labels'])(jj,kk) = a_table.BundleName(kk);
        end
    end
end


%% Plot full results
h = 'right';
close all
cost_fxn = 10;
figure('Position',[100,100,2000,1200]);
for ii = 1:numel(atlas_names)
    a_name = atlas_names{ii};
    d = overlap.(a_name)(:,:,cost_fxn);
    l = overlap.([a_name '_labels'])(1,:);
    subplot(2,3,ii)
    %subtightplot(2,3,ii,[0.05,0.05])
    %bar(d')
    bar(atanh(d'));
    %bar(normalize(atanh(d),2)');
    ax = gca;
    ax.XTick = 1:numel(l);
    ax.XTickLabelRotation = 270;
    ax.XTickLabel = l;
    ax.TickLabelInterpreter = 'none';
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    grid on
    legend(['Digit ROI Projections (' h ')'], ['Letter ROI Projections (' h ')']);
    title(['Track density overlap with ' a_name ' bundles']);
    ax.YLabel.String = 'Correlation between images (Fisher Z)';
    %ax.YLim = [0.5 4.5];
    ax.Title.FontSize = 15;
end

export_fig(['tract_overlap_difference_' h '_all'],'-png','-m2');

%% Plot reduced results
cutoff = 0.15;
cost_fxn = 10;
figure('Position',[100,100,2000,1200]);
for ii = 1:numel(atlas_names)
    a_name = atlas_names{ii};
    d = overlap.(a_name)(:,:,cost_fxn);
    l = overlap.([a_name '_labels'])(1,:);
    
    indkeep = sum(d>=cutoff,1) ~= 0;
    d = d(:,indkeep);
    l = l(:,indkeep);
    
    subplot(2,3,ii)
    %subtightplot(2,3,ii,[0.05,0.05])
    %bar(d')
    bar(atanh(d'));
    %bar(normalize(atanh(d),2)');
    ax = gca;
    ax.XTick = 1:numel(l);
    ax.XTickLabelRotation = 330;
    ax.XTickLabel = l;
    ax.TickLabelInterpreter = 'none';
    ax.FontSize = 13;
    ax.FontWeight = 'bold';
    grid on
    legend(['Digit ROI Projections (' h ')'], ['Letter ROI Projections (' h ')']);
    title(['Track density overlap with ' a_name ' bundles']);
    ax.YLabel.String = 'Correlation between images (Fisher Z)';
    ax.YLim = [0 1.1];
    ax.Title.FontSize = 15;
end

export_fig(['tract_overlap_difference_' h '_reduced'],'-png','-m2');

%% Overlap function
function overlap = calculate_overlap(base,base_ind,source,source_mask)
    base_ind = num2str(base_ind);
    % Compute cost function (registration quality) between base and source using AFNI
    unix(['3dAllineate -quiet -allcostX1D IDENTITY tmp.txt -source_mask ' source_mask ' -base ' base '[' base_ind '] -source ' source]);
    % Open the tmp output file
    fid=fopen('tmp.txt');
    % Set linenum and read in full line
    linenum = 3;
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    C = str2num(cell2mat([C{1,1}]));
    % Record the value of interest (e.g. 13 = lpc+)
    overlap = C;
end

