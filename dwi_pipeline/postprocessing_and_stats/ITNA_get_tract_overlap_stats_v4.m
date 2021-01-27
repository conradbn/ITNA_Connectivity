%% Get overlap of track density images (TDIs) and largescale tracts/bundles
% This script calculates overlap statistics of subject-level TDI maps (in
% MNI space), with all tracts/bundles defined in the (4D) Pandora White
% Matter Atlas @ https://www.nitrc.org/projects/pandora_atlas

purge 

atlas_names = {'Xtract'};% {'AFQ','AFQclipped' 'Recobundles' 'TractSeg' 'Tracula'} 
atlas_dir = '/Users/benconrad/Desktop/Pandora_Atlas';
mask = [atlas_dir '/MNI152_2009_template_mask.nii.gz'];

% Get filenames for each subject
process_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Tract_Overlap/MNI_TDI_maps';
cd(process_dir);
niis = dir('*.MNI.TDI.6mm.nii.gz');

%% Run across subjects/atlases
% % Set thresholds
% tract_thresh = 0.75;
% tdi_thresh = 5;
% % Loop
% for ii = 1:numel(atlas_names)
%     % Read in atlas information
%     a_name = atlas_names{ii};
%     a_table = readtable([atlas_dir '/' a_name '/' a_name '_info.csv']);
%     a_nii = [atlas_dir '/' a_name '/' a_name '.nii'];
%     tracts = niftiread(a_nii); tracts = single(tracts);
%     % Loop through tracts
%     for kk = 1:size(a_table,1)
%         tract = tracts(:,:,:,kk);
%         % Loop through subjects
%         parfor jj = 1:numel(niis)
%             disp(['--------- WORKING ON '...
%                 ' atlas ' num2str(ii) ' of ' num2str(numel(atlas_names))...
%                 ' / tract ' num2str(kk) ' of ' num2str(size(a_table,1))...
%                 ' / subject ' num2str(jj) ' of ' num2str(numel(niis))]);
%             % Load subject TDI map
%             tdi = niftiread(niis(jj).name);
%             % Run function
%             dice(jj) =  calculate_overlap_binary(tract,tdi,tract_thresh,tdi_thresh);
%         end
%         % Add dice values to structure
%         dice_all.(a_name)(kk,:) = dice; clear dice
%     end
%     % Add tract label names to structure
%     dice_all.([a_name '_labels'])(:) = a_table.BundleName';
% end
% save([process_dir '/../dice_Xtract_allsubs.mat']);

%% Plot full results
clear all; close all
process_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Tract_Overlap/MNI_TDI_maps';
cd(process_dir);
load([process_dir '/../dice_Xtract_allsubs.mat']);


%
for ii = 1:numel(atlas_names)
    a_name = atlas_names{ii};
    
    dice = dice_all.(a_name)(:,:)';
    labels = dice_all.([a_name '_labels'])(1,:);
    
    % Recombine the data so it is in a more useful structure
    DigL = dice(1:3:end,:);
    DigR = dice(2:3:end,:);
    LetL = dice(3:3:end,:);
    dice_Combined = cat(3,DigL,LetL,DigR);
    
    % Convert right hemisphere tracts to left, only keep the left values
    for kk = 1:numel(labels)
        if contains(labels{kk},'right')
            rm_label(kk) = true;
            % Get the mirrored DigR data and place in left hemi column
            mirror_label = strrep(labels{kk},'right','left');
            ind = strcmp(labels,mirror_label);
            dice_Combined(:,ind,3) = DigR(:,kk);
        else 
            rm_label(kk) = false;
        end
    end
    
    dice_Combined(:,rm_label,:) = [];
    labels(rm_label) = [];
    
    % Find if significant above zero overlap for any of the ROIs
    cutoff = 0.1;
    for kk = 1:size(dice_Combined,2)
        temp = squeeze(dice_Combined(:,kk,:));
        if any(median(temp) > cutoff)
            ind_keep(kk) = true;
        else
            ind_keep(kk) = false;
        end
    end
    dice_Reduced = dice_Combined(:,ind_keep,:);
    labels_Reduced = labels(1,ind_keep); 
    
    % Reorder based on descending (total) degree of overlap
    %dice_Sum = sum(squeeze(sum(dice_Reduced)),2);
    %dice_Max = max(squeeze(median(dice_Reduced)),[],2);
    dice_Sum_DigL = squeeze(sum(dice_Reduced));
    [~,new_order] = sort(dice_Sum_DigL(:,1),'descend');
    dice_Reduced = dice_Reduced(:,new_order,:);
    labels_Reduced = labels_Reduced(new_order);
    
    %% Run Ttests
    for tt = 1:size(dice_Reduced,2)
        [~,p,ci,t] = ttest(dice_Reduced(:,tt,1),dice_Reduced(:,tt,2));
        stats(tt,1) = t.tstat;
        stats(tt,2) = p;
    end
    stats(:,2) = stats(:,2)./size(dice_Reduced,2); % Bonferoni correction
    
    %% Plot Raw Dice
    close all
    figure('Position',[100,100,2500,700]);
    % Line at zero
    yline(0,'LineWidth',3);
    hold on
    ax = gca;
    % Add patches
    for kk = 1.5:2:numel(labels_Reduced)
        patch([kk kk+1 kk+1 kk],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[.8 .8 .8],'FaceAlpha',0.5,'EdgeColor','none');
    end
    % Custom Boxplots
    import iosr.statistics.*
    bp = boxPlot(dice_Reduced);
    %view([-270 90])
    bp.boxColor = {[.9 .3 .3] [.3 .3 .9] [.9 .7 .7]};
    bp.lineWidth = 3;
    bp.boxWidth = .22;
    bp.lineColor = {[0 0 0] [0 0 0] [0 0 0]};
    bp.medianColor = {[0 0 0] [0 0 0] [0 0 0]};
%     bp.showMean = true;
%     bp.meanMarker = '.';
%     bp.meanColor = 'r';
%     bp.meanSize = 30;
    % Scatter/outliers
    bp.scatterColor = [.4 .4 .4];
    bp.scatterMarker = '.';
    bp.scatterSize = 200;
    bp.showScatter = true;
    bp.symbolColor = [.4 .4 .4];
    bp.symbolMarker = '.';
    bp.showOutliers = true;
    bp.outlierSize = 200;

    grid on
    ax = gca;
    ax.FontSize = 25;
    %ax.YLim = [0.1 0.5];
    ax.GridAlpha = .2;
    ax.GridLineStyle = '--';
    ax.LineWidth = 2;
    % Title
    ax.Title.String = ['Streamline density overlap with ' a_name ' bundles'];
    ax.Title.Interpreter = 'none';
    ax.Title.FontSize = 30;
    % X Ticks
   %ax.XTick = 1:numel(labels_Reduced)*2;
    ax.XTickLabel = labels_Reduced;
    ax.XTickLabelRotation = 330;
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
    low = min(dice_Reduced(:));
    high = max(dice_Reduced(:));
    ax.YLim = [0 high+.05];
    export_fig(['bundle_overlap_' a_name],'-png','-m2','-transparent');

    
    
%% Plot Dice Diff
%     d_Diff = d_Combined(:,:,1) - d_Combined(:,:,2);
%     figure('Position',[100,100,2000,1200]);
%     % Line at zero
%     yline(0,'LineWidth',3);
%     hold on
%     ax = gca;
%     % Add patches
%     for kk = 1.5:2:numel(labels)
%         patch([kk kk+1 kk+1 kk],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[.8 .8 .8],'FaceAlpha',0.5,'EdgeColor','none');
%     end
%     % Custom Boxplots
%     import iosr.statistics.*
%     bp = boxPlot(d_Diff);
%     view([-270 90])
% %     bp.boxColor = {[.9 .3 .3] [.3 .3 .9]};
%     bp.boxColor = {[.6 .3 .6]};
%     bp.lineWidth = 2;
%     bp.boxWidth = .35;
% %     bp.lineColor = {[.9 .3 .3] [.3 .3 .9]};
% %     bp.medianColor = {[.9 .3 .3] [.3 .3 .9]};
%     bp.lineColor = {[.6 .3 .6]};
%     bp.medianColor = {[.6 .3 .6]};
%     bp.showMean = true;
%     bp.meanMarker = '.';
%     bp.meanColor = 'k';
%     bp.meanSize = 30;
%     bp.showOutliers = false;
%     % bp.symbolMarker = '.'
%     % bp.symbolColor = [.8 .8 .8]
%     grid on
%     ax = gca;
%     ax.FontSize = 17;
%     %ax.YLim = [0.1 0.5];
%     ax.GridAlpha = .2;
%     ax.GridLineStyle = '--';
%     % Title
%     ax.Title.String = ['Bundle overlap difference between Digit and Letter ROIs'];
%     ax.Title.Interpreter = 'none';
%     ax.Title.FontSize = 30;
%     % X Ticks
%     ax.XTick = 1:numel(labels);
%     ax.XTickLabel = labels;
%     ax.XAxis.FontSize = 12;
%     ax.XAxis.FontWeight = 'bold';
%     % Axis labels
%     ax.XLabel.String = 'Bundle';
%     ax.XLabel.FontSize = 30;
%     ax.XLabel.FontWeight = 'bold';
%     ax.YLabel.String = 'Difference in Dice coefficient';
%     ax.YLabel.FontSize = 30;
%     ax.YLabel.FontWeight = 'bold';
%     ax.TickLabelInterpreter = 'none';
%     
%     % Reset Yaxis 
%     low = min(d_Diff(:));
%     high = max(d_Diff(:));
%     ax.YLim = [low-.05 high+.05];
%     
%     export_fig(['bundle_overlap_differences_' h '_' a_name],'-png','-m2');
%     
%     
% %     % TTests
% %     [htest,p,ci,stat] = ttest(d_Digit,d_Letter,'Alpha',0.05/size(d,2)); % Bonferoni corrected
% %     ind = 1;
% %     for kk = 1:numel(htest)
% %         if htest(kk) == 1
% %             t_sig(ind) = stat.tstat(kk);
% %             l2(ind) = l(kk);
% %             ind = ind + 1;
% %         else
% %             continue
% %         end
% %     end
% %     figure('Position',[100,100,2000,1200]);
% %     hold on
% %     for kk = 1:numel(t_sig)
% %         if t_sig(kk) > 0
% %             b = bar(kk,t_sig(kk));
% %             b.FaceColor = [.9 .3 .3];
% %         else
% %             b = bar(kk,t_sig(kk));
% %             b.FaceColor = [.3 .3 .9];
% %         end
% %     end
% %     
% %     ax = gca;
% %     grid on
% %     ax = gca;
% %     ax.FontSize = 17;
% %     %ax.YLim = [0.1 0.5];
% %     ax.GridAlpha = .2;
% %     ax.GridLineStyle = '--';
% %     % Title
% %     ax.Title.String = [a_name ' bundles showing significant difference in overlap with ' h ' hemisphere ROI projections'];
% %     ax.Title.Interpreter = 'none';
% %     ax.Title.FontSize = 30;
% %     % X Ticks
% %     ax.XTick = 1:numel(l2);
% %     ax.XTickLabelRotation = 270;
% %     ax.XTickLabel = l2;
% %     ax.XTickLabelRotation = 270;
% %     ax.XAxis.FontSize = 12;
% %     ax.XAxis.FontWeight = 'bold';
% %     % Axis labels
% %     ax.XLabel.String = 'Bundle';
% %     ax.XLabel.FontSize = 30;
% %     ax.XLabel.FontWeight = 'bold';
% %     ax.YLabel.String = 'T-statistic of Difference in overlap (paired test)';
% %     ax.YLabel.FontSize = 30;
% %     ax.YLabel.FontWeight = 'bold';
% %     ax.TickLabelInterpreter = 'none';
% %     
% %     %export_fig(['bundle_overlap_differences_' h '_' a_name],'-png','-m2');
% %     
% %     clear t_sig l2 l b
%     
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
    


