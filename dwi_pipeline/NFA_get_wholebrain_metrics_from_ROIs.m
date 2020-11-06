%% Extract length and density metrics within ROIs, from subject-level wholebrain maps
% Plot and run statistical comparisons

purge;
cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface')
data_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Data';
subs = dir('1*');

roi_l = {'std.141.lh.PollackandPrice19_Dp-Da.MNI152.votc.inflated.14mm_diam.1.1D'    
         'std.141.lh.PollackandPrice19_Lp-La.MNI152.votc.inflated.14mm_diam.1.1D'};
roi_r = {'std.141.rh.PollackandPrice19_Dp-Da.MNI152.votc.inflated.14mm_diam.1.1D'
         'std.141.rh.PollackandPrice19_Lp-La.MNI152.votc.inflated.14mm_diam.1.1D'};
data_l = {'all_subs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.log.niml.dset'
          'all_subs_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset'};
data_r = {'all_subs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.log.niml.dset'
          'all_subs_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset'};


for ii = 1:numel(subs)
    sub = subs(ii).name;
    suma_dir = ['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sub '_proc/' sub '.freesurfer/SUMA'];
    for jj = 1:numel(data_l)
        % Load data for this subject
        dall = afni_niml_readsimple([data_dir '/' data_l{jj}]);
        dall = dall.data(:,ii);
        for kk = 1:numel(roi_l)
            % Load the subject's circular ROI
            m = readmatrix([suma_dir '/' roi_l{kk}],'NumHeaderLines',2,'FileType','text');
            % Convert indices
            inds = m(:,1)+1;
            d = dall(inds);
            out_l_mean(ii,kk,jj) = mean(d);
            out_l_std(ii,kk,jj) = std(d);
        end
    end
    for jj = 1:numel(data_r)
        % Load data for this subject
        dall = afni_niml_readsimple([data_dir '/' data_r{jj}]);
        dall = dall.data(:,ii);
        for kk = 1:numel(roi_r)
            % Load the subject's circular ROI
            m = readmatrix([suma_dir '/' roi_r{kk}],'NumHeaderLines',2,'FileType','text');
            % Convert indices
            inds = m(:,1)+1;
            d = dall(inds);
            out_r_mean(ii,kk,jj) = mean(d);
            out_r_std(ii,kk,jj) = std(d);
        end
    end
end

%% Convert to table for easier later parsing
T = table();
T.("Dp-Lp_seed_TDIends_lh") = out_l_mean(:,1,1);
T.("Dp-Lp_seed_Length_lh") = out_l_mean(:,1,2);
T.("Dp-Lp_seed_TDIends_rh") = out_r_mean(:,2,1);
T.("Dp-Lp_seed_Length_rh") = out_r_mean(:,2,2);
save('TDI_endpoint_density_and_length.mat','T');

%% Ttests

[~,p,~,t] = ttest(out_l_mean(:,1,1),out_l_mean(:,2,1));
disp(['Density LH - D v L - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);
[~,p,~,t] = ttest(out_l_mean(:,1,2),out_l_mean(:,2,2));
disp(['Length LH - D v L - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);
[~,p,~,t] = ttest(out_r_mean(:,1,1),out_r_mean(:,2,1));
disp(['Density RH - D v L - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);
[~,p,~,t] = ttest(out_r_mean(:,1,2),out_r_mean(:,2,2));
disp(['Length RH - D v L - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);


%% Plot
close all
length = cat(3,out_l_mean(:,:,2),out_r_mean(:,:,2));
length = permute(length,[1,3,2]);
make_boxplot(length)
ax = gca;
ax.Title.String = 'Streamline Length - Digit versus Letter ROIs';
ax.YLabel.String = 'Mean, Weighted Length (mm)';
ax.XTickLabel = {'Left Hemisphere' 'Right Hemisphere'};
ax.FontSize = 20;
ax.GridAlpha = 0.8;
ax.GridLineStyle = ':';
export_fig('Final_Figs/wholebrain_length','-png','-m2','-transparent');

density = cat(3,out_l_mean(:,:,1),out_r_mean(:,:,1));
density = permute(density,[1,3,2]);
make_boxplot(density)
ax = gca;
ax.Title.String;
ax.Title.String = 'Connection Density - Digit versus Letter ROIs';
ax.YLabel.String = 'Mean, Weighted Density - log_1_0(mm^2)';
ax.XTickLabel = {'Left Hemisphere' 'Right Hemisphere'};
ax.FontSize = 20;
ax.GridAlpha = 0.8;
ax.GridLineStyle = ':';
export_fig('Final_Figs/wholebrain_density','-png','-m2','-transparent');

function make_boxplot(d)
    import iosr.statistics.*
    figure('Position',[100 100 900 600])
    % Boxplot
    bp = boxPlot(d);
    bp.boxColor = {[.9 .3 .3] [.3 .3 .9]};
    bp.lineWidth = 2;
    bp.boxWidth = .35;
    grid on;
    % Mean/median
    bp.showMean = true;
    bp.meanMarker = '.';
    bp.meanColor = 'k';
    bp.meanSize = 30;
    bp.lineColor = {[.9 .3 .3] [.3 .3 .9]};
    bp.medianColor = {[.9 .3 .3] [.3 .3 .9]};

    bp.scatterColor = 'k';
    bp.scatterMarker = 'o';
    bp.scatterSize = 30;
    bp.showScatter = true;
    bp.symbolColor = 'k';
    bp.symbolMarker = 'o';
    bp.showOutliers = true;
end


