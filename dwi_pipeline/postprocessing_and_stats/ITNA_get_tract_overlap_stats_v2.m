%% Get overlap of TDI and largescale tracts
% This script calculates overlap statistics of subject-level TDI maps (in
% MNI space), with all tracts/bundles defined in the (4D) Pandora White
% Matter Atlas @ https://www.nitrc.org/projects/pandora_atlas

purge 

atlas_names = {'AFQ'};% 'Xtract'};%'AFQclipped' 'Recobundles' 'TractSeg' 'Tracula' 
atlas_dir = '/Volumes/NBL_Projects/Price_NFA/NFA_DWI/Pandora_Atlas';%'/Users/benconrad/Downloads/NITRC-multi-file-downloads'; 
mask = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152/MNI152_2009_template_mask.nii.gz';

% Get filenames for each subject
%cd /Users/nbl_imac2/Desktop/Price_NFA_Tractography_MNI152/
cd /Volumes/BensHD_2020/Price_NFA_Tractography_MNI152
%niis = subdir('Group_Average_*p-*a_rh_fwhm4.nii.gz'); %'tracks_ss3t_*_combined_Dp-Da_MNI_TDI_fwhm4.nii.gz');
%niis = subdir('Digit-Letter_TDI_rh_4mm_3dttest++_Coef_masked_*.nii'); 
niis = subdir('tracks_ss3t_*_combined_*p-*a.rh.MNI.TDI.smooth6mm.nii.gz');

%% Run across subjects/atlases
for ii = 1:numel(atlas_names)
    a_name = atlas_names{ii};
    a_table = readtable([atlas_dir '/' a_name '/' a_name '_info.csv']);
    a_nii = [atlas_dir '/' a_name '/' a_name '.nii.gz']; %load_untouch_nii_gz([atlas_dir '/' a_name '/' a_name '.nii.gz']);
    for jj = 1:numel(niis)
        n = niis(jj).name;
        parfor kk = 1:size(a_table,1)
%             disp(['--------- WORKING ON '...
%                   ' atlas ' num2str(ii) ' of ' num2str(numel(atlas_names)) ...
%                   ' / subject ' num2str(jj) ' of ' num2str(numel(niis))...
%                   ' / tract ' num2str(kk) ' of ' num2str(size(a_table,1))])
            % Setup inputs
            base = a_nii;
            base_ind = kk-1;
            source = n;
            source_mask = mask;
            %cost_fxn = 13; % LPC+, combination metric
            % Run function
            ov(kk,:) = calculate_overlap(base,base_ind,source,source_mask);
            ovname(kk) = a_table.BundleName(kk);
        end
       overlap.(a_name)(jj) = ov; clear ov
       overlap.([a_name '_labels'])(jj) = ovname;    
    end
end


%% Plot full results
cd /Volumes/NBL_Projects/Price_NFA/NFA_DWI/Group_StatisticalTests
load('/Volumes/NBL_Projects/Price_NFA/NFA_DWI/Group_StatisticalTests/overlap_AFQ_Xtract_allsubs_lh.mat');
h = 'left';
close all
cost_fxn = 14; %10 =Pearson correlation
for ii = 1:numel(atlas_names)
    a_name = atlas_names{ii};
    d = overlap.(a_name)(:,:,cost_fxn);
    l = overlap.([a_name '_labels'])(1,:);
    d_Digit = d(1:2:end,:);
    d_Letter = d(2:2:end,:);
    d_Combined = cat(3,d_Digit,d_Letter);
    % Plot
    figure('Position',[100,100,2000,1200]);
    % Line at zero
    yline(0,'LineWidth',3);
    hold on
    ax = gca;
    % Add patches
    for kk = 1.5:2:numel(l)
        patch([kk kk+1 kk+1 kk],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[.8 .8 .8],'FaceAlpha',0.5,'EdgeColor','none');
    end
    
    
    % Custom Boxplots
    import iosr.statistics.*
    bp = boxPlot(d_Combined);
    bp.boxColor = {[.9 .3 .3] [.3 .3 .9]};
    bp.lineWidth = 2;
    bp.boxWidth = .35;
    bp.lineColor = {[.9 .3 .3] [.3 .3 .9]};
    bp.medianColor = {[.9 .3 .3] [.3 .3 .9]};
    bp.showMean = true;
    bp.meanMarker = '.';
    bp.meanColor = 'k';
    bp.meanSize = 30;
    bp.showOutliers = false;
    % bp.symbolMarker = '.'
    % bp.symbolColor = [.8 .8 .8]
    grid on
    ax = gca;
    ax.FontSize = 17;
    %ax.YLim = [0.1 0.5];
    ax.GridAlpha = .2;
    ax.GridLineStyle = '--';
    % Title
    ax.Title.String = ['Streamline density overlap with ' a_name ' bundles for ' h ' hemisphere ROIs'];
    ax.Title.Interpreter = 'none';
    ax.Title.FontSize = 30;
    % X Ticks
    ax.XTick = 1:numel(l)*2;
    ax.XTickLabelRotation = 270;
    ax.XTickLabel = l;
    ax.XTickLabelRotation = 270;
    ax.XAxis.FontSize = 12;
    ax.XAxis.FontWeight = 'bold';
    % Axis labels
    ax.XLabel.String = 'Bundle';
    ax.XLabel.FontSize = 30;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.String = 'Correlation between images (Fisher Z)';
    ax.YLabel.FontSize = 30;
    ax.YLabel.FontWeight = 'bold';
    
    % Reset Yaxis 
    low = min(d_Combined(:));
    high = max(d_Combined(:));
    ax.YLim = [low-.05 high+.05];
    %export_fig(['bundle_overlap_' h '_' a_name],'-png','-m2');
    
    % TTests
    [htest,p,ci,stat] = ttest(d_Digit,d_Letter,'Alpha',0.05/size(d,2)); % Bonferoni corrected
    
    ind = 1;
    for kk = 1:numel(htest)
        if htest(kk) == 1
            t_sig(ind) = stat.tstat(kk);
            l2(ind) = l(kk);
            ind = ind + 1;
        else
            continue
        end
    end
    figure('Position',[100,100,2000,1200]);
    hold on
    for kk = 1:numel(t_sig)
        if t_sig(kk) > 0
            b = bar(kk,t_sig(kk));
            b.FaceColor = [.9 .3 .3];
        else
            b = bar(kk,t_sig(kk));
            b.FaceColor = [.3 .3 .9];
        end
    end
    
    ax = gca;
    grid on
    ax = gca;
    ax.FontSize = 17;
    %ax.YLim = [0.1 0.5];
    ax.GridAlpha = .2;
    ax.GridLineStyle = '--';
    % Title
    ax.Title.String = [a_name ' bundles showing significant difference in overlap with ' h ' hemisphere ROI projections'];
    ax.Title.Interpreter = 'none';
    ax.Title.FontSize = 30;
    % X Ticks
    ax.XTick = 1:numel(l2);
    ax.XTickLabelRotation = 270;
    ax.XTickLabel = l2;
    ax.XTickLabelRotation = 270;
    ax.XAxis.FontSize = 12;
    ax.XAxis.FontWeight = 'bold';
    % Axis labels
    ax.XLabel.String = 'Bundle';
    ax.XLabel.FontSize = 30;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.String = 'T-statistic of Difference in overlap (paired test)';
    ax.YLabel.FontSize = 30;
    ax.YLabel.FontWeight = 'bold';
    
    %export_fig(['bundle_overlap_differences_' h '_' a_name],'-png','-m2');
    
    clear t_sig l2 l b
    
end




%     %subplot(2,3,ii)
%     %subtightplot(2,3,ii,[0.05,0.05])
%     %bar(d')
%     notBoxPlot(atanh(d_Reorder));
%     
%     %bar(normalize(atanh(d),2)');
%     ax = gca;
%     ax.XTick = 1:numel(l);
%     ax.XTickLabelRotation = 270;
%     ax.XTickLabel = l;
%     ax.TickLabelInterpreter = 'none';
%     ax.FontSize = 10;
%     ax.FontWeight = 'bold';
%     grid on
%     legend(['Digit ROI Projections (' h ')'], ['Letter ROI Projections (' h ')']);
%     title(['Track density overlap with ' a_name ' bundles']);
%     ax.YLabel.String = 'Correlation between images (Fisher Z)';
%     %ax.YLim = [0.5 4.5];
%     ax.Title.FontSize = 15;
%     
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

