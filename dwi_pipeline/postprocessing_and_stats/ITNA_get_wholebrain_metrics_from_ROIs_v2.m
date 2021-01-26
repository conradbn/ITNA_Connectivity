%% Extract length and density metrics within ROIs, from subject-level wholebrain maps
% Plot and run statistical comparisons
% NOTE - This code has been updated to use files on the NBL_Projects server
% space, and focuses on only the Digit L, Letter L, and Digit R (math) ROIs

purge;
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData')
subs = dir('1*');
data_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/AllSubs_dsets';


roi = {'std.141.lh.PollackandPrice19_Dp-Da.MNI152.votc.inflated.14mm_diam.1.1D'    
       'std.141.lh.PollackandPrice19_Lp-La.MNI152.votc.inflated.14mm_diam.1.1D'
       'std.141.rh.PP19_Dp-Da_math.MNI152.votc.inflated.14mm_diam.1.1D'};

% data = {'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.lh.6mm.log.niml.dset'
%         'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.lh.6mm.niml.dset'
%         'AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.rh.6mm.log.niml.dset'
%         'AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.rh.6mm.niml.dset'};
 

for kk = 1:numel(roi)
    if contains(roi{kk},'lh')
        h = 'lh';
    else 
        h = 'rh';
    end
    % Load full catenated dataset 
    dens = afni_niml_readsimple([data_dir '/AllSubs_tracks_ss3t_50M.wholebrain_TDI_ends.norm.al2anat.' h '.6mm.niml.dset']);
    leng = afni_niml_readsimple([data_dir '/AllSubs_tracks_ss3t_50M.wholebrain_length_map.al2anat.' h '.6mm.niml.dset']);
    for ii = 1:numel(subs)
        sub = subs(ii).name;
        sub = strrep(sub,'_proc','');
        disp(sub);
        suma_dir = ['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sub '_proc/' sub '.freesurfer/SUMA'];
        % Load the subject's circular ROI
        m = readmatrix([suma_dir '/' roi{kk}],'NumHeaderLines',2,'FileType','text');
        d = dens.data(:,ii);
        l = leng.data(:,ii);
        % Convert indices
        inds = m(:,1)+1;
        dinds = d(inds);
        linds = l(inds);
        out_mean(kk,ii,1) = sum(dinds);
        out_mean(kk,ii,2) = mean(linds);
        out_std(kk,ii,1) = std(dinds);
        out_std(kk,ii,1) = std(linds);
    end
end

save(['/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Density_Length_Comparisons/'...
     'density_length_comparison_initial_workspace.mat']);


%% Convert to table, run ttests, and plot
clear all
out_dir = '/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Density_Length_Comparisons';
load([out_dir '/density_length_comparison_initial_workspace.mat']);

T = table();
T.("DigLH_density") = log10(out_mean(1,:,1)');
T.("DigLH_length") = out_mean(1,:,2)';
T.("LetLH_density") = log10(out_mean(2,:,1)');
T.("LetLH_length") = out_mean(2,:,2)';
T.("DigRH_density") = log10(out_mean(3,:,1)');
T.("DigRH_length") = out_mean(3,:,2)';


% Ttests
[~,p,~,t] = ttest(T.DigLH_density,T.LetLH_density);
disp(['Density - DigLH v LetLH - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);
[~,p,~,t] = ttest(T.DigLH_density,T.DigRH_density);
disp(['Density - DigLH v DigRH - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);
[~,p,~,t] = ttest(T.DigLH_length,T.LetLH_length);
disp(['Length - DigLH v LetLH - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);
[~,p,~,t] = ttest(T.DigLH_length,T.DigRH_length);
disp(['Length - DigLH v DigRH - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);

[~,p,~,t] = ttest(T.DigRH_density,T.LetLH_density);
disp(['Density - DigRH v LetLH - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);
[~,p,~,t] = ttest(T.DigRH_length,T.LetLH_length);
disp(['Length - DigRH v LetLH - tstat = ' num2str(t.tstat) ', p = ' num2str(p)]);

 
% Plot
close all
length = [T.DigLH_length,T.LetLH_length,T.DigRH_length];
make_boxplot(length)
ax = gca;
ax.Title.String = 'Streamline Length';
ax.YLabel.String = 'Mean, Weighted Length (mm)';
ax.XTickLabel = {'Digit L' 'Letter L' 'Digit R'};
ax.FontSize = 30;
ax.GridAlpha = 0.8;
ax.GridLineStyle = ':';
ax.Title.FontSize = 40;
ax.YLabel.FontSize = 32;
export_fig([out_dir '/wholebrain_length'],'-png','-m2','-transparent');

density = [T.DigLH_density,T.LetLH_density,T.DigRH_density];
make_boxplot(density)
ax = gca;
ax.Title.String;
ax.Title.String = 'Fiber Density';
ax.YLabel.String = 'Total cross-sectional area - log_1_0(mm^2)';
ax.XTickLabel = {'Digit L' 'Letter L' 'Digit R'};
ax.FontSize = 30;
ax.GridAlpha = 0.8;
ax.GridLineStyle = ':';
ax.Title.FontSize = 40;
ax.YLabel.FontSize = 31;

ax.YLim = [-.2 .7];
export_fig([out_dir '/wholebrain_density'],'-png','-m2','-transparent');

function make_boxplot(d)
    import iosr.statistics.*
    figure('Position',[100 100 700 550])
    % Boxplot
    bp = boxPlot(d);
    bp.boxColor = {[.8 .8 .8]};
    bp.lineWidth = 2;
    bp.boxWidth = .75;
    grid on;
    % Mean/median
    bp.showMean = true;
    bp.meanMarker = '.';
    bp.meanColor = [.9 .3 .3];
    bp.meanSize = 50;
    bp.lineColor = {[.9 .3 .3]};
    bp.medianColor = {[.9 .3 .3]};

    bp.scatterColor = 'k';
    bp.scatterMarker = '.';
    bp.scatterSize = 400;
    bp.showScatter = true;
    bp.symbolColor = 'k';
    bp.symbolMarker = '.';
    bp.outlierSize = 400;
    bp.showOutliers = true;
end


