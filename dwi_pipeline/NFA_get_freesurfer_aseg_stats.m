
purge;
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData'); 
sub_dirs = dir('*_proc');

%% Create string containing all file paths/names
file_list = [];
for ii = 1:numel(sub_dirs)
    sub = strsplit(sub_dirs(ii).name,'_');
    sub = sub{1};
    s{ii} = sub;
    % Stats file path
    stats_file = ['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' sub '_proc/' sub '.freesurfer/stats/aseg.stats'];
    file_list = [file_list ' ' stats_file];
end

%% Make stats table file
cd('/Volumes/BensHD_2020/Price_NFA_Tractography_MNI152');
unix(['asegstats2table --inputs ' file_list ' --meas volume --tablefile aseg_stats.txt'])

%% Read table and add subject ID 
T = readtable('aseg_stats.txt');
T.subID = s';

%% Write to text
writetable(T,'aseg_stats_subid.csv')