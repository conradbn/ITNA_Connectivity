%% Get overlap of track density images (TDIs) and largescale tracts/bundles
% This script calculates overlap statistics of subject-level TDI maps (in
% MNI space), with all tracts/bundles defined in the (4D) Pandora White
% Matter Atlas @ https://www.nitrc.org/projects/pandora_atlas

purge 

atlas_names = {'AFQ','Xtract'};% 'Xtract'};%'AFQclipped' 'Recobundles' 'TractSeg' 'Tracula' 
atlas_dir = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152/Pandora_Atlas';%'/Users/benconrad/Downloads/NITRC-multi-file-downloads'; 
mask = '/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152/MNI152_2009_template_mask.nii.gz';

% Get filenames for each subject
cd /Volumes/BensHD_2020/Price_NFA_Tractography_MNI152
niis = subdir('tracks_ss3t_*_combined_*p-*a.lh.MNI.TDI.smooth6mm.nii.gz');

%% Run across subjects/atlases
% Set thresholds
tract_thresh = 0.25;
tdi_thresh = 5;
% Loop
for ii = 1:numel(atlas_names)
    % Read in atlas information
    a_name = atlas_names{ii};
    a_table = readtable([atlas_dir '/' a_name '/' a_name '_info.csv']);
    a_nii = [atlas_dir '/' a_name '/' a_name '.nii'];
    tracts = niftiread(a_nii); tracts = single(tracts);
    % Loop through tracts
    for kk = 1:size(a_table,1)
        tract = tracts(:,:,:,kk);
        % Loop through subjects
        parfor jj = 1:numel(niis)
            disp(['--------- WORKING ON '...
                ' atlas ' num2str(ii) ' of ' num2str(numel(atlas_names))...
                ' / tract ' num2str(kk) ' of ' num2str(size(a_table,1))...
                ' / subject ' num2str(jj) ' of ' num2str(numel(niis))]);
            % Load subject TDI map
            tdi = niftiread(niis(jj).name);
            % Run function
            dice(jj) =  calculate_overlap_binary(tract,tdi,tract_thresh,tdi_thresh);
        end
        % Add dice values to structure
        dice_all.(a_name)(kk,:) = dice; clear dice
    end
    % Add tract label names to structure
    dice_all.([a_name '_labels'])(:) = a_table.FileName';
end


%% Plot full results
%cd /Volumes/NBL_Projects/Price_NFA/NFA_DWI/Group_StatisticalTests
%load('/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152/dice_AFQ_Xtract_allsubs_lh.mat');
save('/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152/dice_AFQ_Xtract_allsubs_lh.mat');
h = 'left';
close all

for ii = 1:numel(atlas_names)
    a_name = atlas_names{ii};
    d = dice_all.(a_name)(:,:)';
    l = dice_all.([a_name '_labels'])(1,:);
    d_Digit = d(1:2:end,:);
    d_Letter = d(2:2:end,:);
    d_Combined = cat(3,d_Digit,d_Letter);
    
    % Plot Raw Dice
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
    ax.YLabel.String = 'Dice coefficient between thresholded/binarized images';
    ax.YLabel.FontSize = 30;
    ax.YLabel.FontWeight = 'bold';
    ax.TickLabelInterpreter = 'none';
    
    % Reset Yaxis 
    low = min(d_Combined(:));
    high = max(d_Combined(:));
    ax.YLim = [low-.05 high+.05];
    export_fig(['bundle_overlap_' h '_' a_name],'-png','-m2');

    
    
    % Plot Dice Diff
    d_Diff = d_Digit - d_Letter;
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
    bp = boxPlot(d_Diff);
%     bp.boxColor = {[.9 .3 .3] [.3 .3 .9]};
    bp.boxColor = {[.6 .3 .6]};
    bp.lineWidth = 2;
    bp.boxWidth = .35;
%     bp.lineColor = {[.9 .3 .3] [.3 .3 .9]};
%     bp.medianColor = {[.9 .3 .3] [.3 .3 .9]};
    bp.lineColor = {[.6 .3 .6]};
    bp.medianColor = {[.6 .3 .6]};
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
    ax.Title.String = ['Bundle overlap difference between Digit and Letter regions'];
    ax.Title.Interpreter = 'none';
    ax.Title.FontSize = 30;
    % X Ticks
    ax.XTick = 1:numel(l);
    ax.XTickLabelRotation = 270;
    ax.XTickLabel = l;
    ax.XTickLabelRotation = 270;
    ax.XAxis.FontSize = 12;
    ax.XAxis.FontWeight = 'bold';
    % Axis labels
    ax.XLabel.String = 'Bundle';
    ax.XLabel.FontSize = 30;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.String = 'Difference in Dice coefficient';
    ax.YLabel.FontSize = 30;
    ax.YLabel.FontWeight = 'bold';
    ax.TickLabelInterpreter = 'none';
    
    % Reset Yaxis 
    low = min(d_Diff(:));
    high = max(d_Diff(:));
    ax.YLim = [low-.05 high+.05];
    
    export_fig(['bundle_overlap_differences_' h '_' a_name],'-png','-m2');
%     % TTests
%     [htest,p,ci,stat] = ttest(d_Digit,d_Letter,'Alpha',0.05/size(d,2)); % Bonferoni corrected
%     ind = 1;
%     for kk = 1:numel(htest)
%         if htest(kk) == 1
%             t_sig(ind) = stat.tstat(kk);
%             l2(ind) = l(kk);
%             ind = ind + 1;
%         else
%             continue
%         end
%     end
%     figure('Position',[100,100,2000,1200]);
%     hold on
%     for kk = 1:numel(t_sig)
%         if t_sig(kk) > 0
%             b = bar(kk,t_sig(kk));
%             b.FaceColor = [.9 .3 .3];
%         else
%             b = bar(kk,t_sig(kk));
%             b.FaceColor = [.3 .3 .9];
%         end
%     end
%     
%     ax = gca;
%     grid on
%     ax = gca;
%     ax.FontSize = 17;
%     %ax.YLim = [0.1 0.5];
%     ax.GridAlpha = .2;
%     ax.GridLineStyle = '--';
%     % Title
%     ax.Title.String = [a_name ' bundles showing significant difference in overlap with ' h ' hemisphere ROI projections'];
%     ax.Title.Interpreter = 'none';
%     ax.Title.FontSize = 30;
%     % X Ticks
%     ax.XTick = 1:numel(l2);
%     ax.XTickLabelRotation = 270;
%     ax.XTickLabel = l2;
%     ax.XTickLabelRotation = 270;
%     ax.XAxis.FontSize = 12;
%     ax.XAxis.FontWeight = 'bold';
%     % Axis labels
%     ax.XLabel.String = 'Bundle';
%     ax.XLabel.FontSize = 30;
%     ax.XLabel.FontWeight = 'bold';
%     ax.YLabel.String = 'T-statistic of Difference in overlap (paired test)';
%     ax.YLabel.FontSize = 30;
%     ax.YLabel.FontWeight = 'bold';
%     ax.TickLabelInterpreter = 'none';
%     
%     %export_fig(['bundle_overlap_differences_' h '_' a_name],'-png','-m2');
%     
%     clear t_sig l2 l b
    
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
% cutoff = 0.15;
% cost_fxn = 10;
% figure('Position',[100,100,2000,1200]);
% for ii = 1:numel(atlas_names)
%     a_name = atlas_names{ii};
%     d = overlap.(a_name)(:,:,cost_fxn);
%     l = overlap.([a_name '_labels'])(1,:);
%     
%     indkeep = sum(d>=cutoff,1) ~= 0;
%     d = d(:,indkeep);
%     l = l(:,indkeep);
%     
%     subplot(2,3,ii)
%     %subtightplot(2,3,ii,[0.05,0.05])
%     %bar(d')
%     bar(atanh(d'));
%     %bar(normalize(atanh(d),2)');
%     ax = gca;
%     ax.XTick = 1:numel(l);
%     ax.XTickLabelRotation = 330;
%     ax.XTickLabel = l;
%     ax.TickLabelInterpreter = 'none';
%     ax.FontSize = 13;
%     ax.FontWeight = 'bold';
%     grid on
%     legend(['Digit ROI Projections (' h ')'], ['Letter ROI Projections (' h ')']);
%     title(['Track density overlap with ' a_name ' bundles']);
%     ax.YLabel.String = 'Correlation between images (Fisher Z)';
%     ax.YLim = [0 1.1];
%     ax.Title.FontSize = 15;
% end
% 
% export_fig(['tract_overlap_difference_' h '_reduced'],'-png','-m2');

%% Dice overlap function
function dice = calculate_overlap_binary(tract,tdi,tract_thresh,tdi_thresh)
    % Load the nifti images
    sthresh = tdi >= tdi_thresh; 
    bthresh = tract >= tract_thresh;
    Seg1 = sthresh(:);
    Seg2 = bthresh(:);
    VoxelsNumber1=sum(Seg1); 
    VoxelsNumber2=sum(Seg2);
    CommonArea=sum(Seg1 & Seg2); 
    dice=(2*CommonArea)/(VoxelsNumber1+VoxelsNumber2);
end
    


